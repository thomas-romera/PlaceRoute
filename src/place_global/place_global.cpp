#include "place_global.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>

#include "density_legalizer.hpp"
#include "net_model.hpp"

namespace coloquinte {

namespace {
/**
 * Mix two placements with a weight: 0 for the first, 1 for the second, or
 * in-between
 */
std::vector<float> blendPlacement(const std::vector<float> &v1,
                                  const std::vector<float> &v2,
                                  float blending) {
  if (blending == 0.0f) {
    return v1;
  }
  if (blending == 1.0f) {
    return v2;
  }
  std::vector<float> ret;
  assert(v1.size() == v2.size());
  ret.reserve(v1.size());
  for (size_t i = 0; i < v1.size(); ++i) {
    ret.push_back((1.0f - blending) * v1[i] + blending * v2[i]);
  }
  return ret;
}
}  // namespace

void GlobalPlacer::place(Circuit &circuit, const ColoquinteParameters &params,
                         const std::optional<PlacementCallback> &callback) {
  params.check();
  std::cout << "Global placement starting" << std::endl;
  auto startTime = std::chrono::steady_clock::now();
  GlobalPlacer pl(circuit, params);
  pl.callback_ = callback;
  pl.run();
  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<float> duration = endTime - startTime;
  std::cout << std::fixed << std::setprecision(2) << "Global placement done in "
            << duration.count() << "s" << std::endl;
  pl.exportPlacement(circuit);
}

GlobalPlacer::GlobalPlacer(Circuit &circuit, const ColoquinteParameters &params)
    : leg_(DensityLegalizer::fromIspdCircuit(
          circuit, params.global.roughLegalization.binSize,
          params.global.roughLegalization.sideMargin)),
      xtopo_(NetModel::xTopology(circuit)),
      ytopo_(NetModel::yTopology(circuit)),
      params_(params),
      circuit_(circuit) {
  circuit_.hasCellSizeUpdate_ = false;
  circuit_.hasNetUpdate_ = false;
  rgen_.seed(params_.seed);
  averageCellLength_ = computeAverageCellSize();
  perCellPenalty_ = computePerCellPenalty();
  auto rlp = params.global.roughLegalization;
  LegalizationModel m = rlp.costModel;
  DensityLegalizer::Parameters legParams;
  legParams.nbSteps = rlp.nbSteps;
  legParams.costModel = rlp.costModel;
  legParams.lineReoptSize = rlp.lineReoptSize;
  legParams.lineReoptOverlap = rlp.lineReoptOverlap;
  legParams.diagReoptSize = rlp.diagReoptSize;
  legParams.diagReoptOverlap = rlp.diagReoptOverlap;
  legParams.squareReoptSize = rlp.squareReoptSize;
  legParams.squareReoptOverlap = rlp.squareReoptOverlap;
  legParams.unidimensionalTransport =
      rlp.unidimensionalTransport && m == LegalizationModel::L1;
  legParams.coarseningLimit = rlp.coarseningLimit;
  if (m == LegalizationModel::L1 || m == LegalizationModel::L2 ||
      m == LegalizationModel::LInf) {
    float dist = leg_.placementArea().width() + leg_.placementArea().height();
    legParams.quadraticPenaltyFactor = rlp.quadraticPenalty / dist;
  }
  leg_.setParams(legParams);
  initDensification();
}

void GlobalPlacer::initDensification() {
  // Capture the true widths so every step's inflation is computed from the
  // baseline (never compounding), and the real sizes can be restored at export.
  baseCellWidth_ = circuit_.cellWidth();
  densificationActive_ = false;
  long long movableArea = 0;
  for (int i = 0; i < circuit_.nbCells(); ++i) {
    if (!circuit_.isFixed(i)) {
      movableArea += circuit_.area(i);
    }
  }
  long long capacity = leg_.totalCapacity();
  baseDensity_ = capacity > 0 ? (double)movableArea / (double)capacity : 1.0;
}

void GlobalPlacer::exportPlacement(Circuit &circuit) const {
  assert(xtopo_.nbCells() == circuit.nbCells());
  assert(ytopo_.nbCells() == circuit.nbCells());
  assert(leg_.nbCells() == circuit.nbCells());
  float w = params_.global.exportBlending;
  std::vector<float> xplace = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yplace = blendPlacement(yPlacementLB_, yPlacementUB_, w);
  exportPlacement(circuit, xplace, yplace);
}

void GlobalPlacer::exportPlacement(Circuit &circuit,
                                   const std::vector<float> &xplace,
                                   const std::vector<float> &yplace) {
  assert((int)xplace.size() == circuit.nbCells());
  assert((int)yplace.size() == circuit.nbCells());

  for (int i = 0; i < circuit.nbCells(); ++i) {
    if (circuit.isFixed(i)) {
      continue;
    }
    circuit.cellX_[i] = std::round(xplace[i] - 0.5 * circuit.placedWidth(i));
    circuit.cellY_[i] = std::round(yplace[i] - 0.5 * circuit.placedHeight(i));
  }
}

float GlobalPlacer::computeAverageCellSize() const {
  float totalDemand = leg_.totalDemand();
  float avgDemand = totalDemand == 0.0 ? 0.0 : totalDemand / leg_.nbCells();
  return std::sqrt(avgDemand);
}

std::vector<float> GlobalPlacer::computePerCellPenalty() const {
  int nbCells = leg_.nbNonEmptyCells();
  float meanArea = leg_.totalDemand() / std::max(1, nbCells);
  std::vector<float> ret;
  ret.reserve(leg_.nbCells());

  for (int i = 0; i < leg_.nbCells(); ++i) {
    ret.push_back(std::pow(leg_.cellDemand(i) / meanArea,
                           params_.global.penalty.areaExponent));
  }
  return ret;
}

std::vector<float> GlobalPlacer::computeIterationPerCellPenalty() {
  std::vector<float> penalty = perCellPenalty_;
  for (float &s : penalty) {
    s *= penalty_;
  }
  if (params_.global.noise > 0.0) {
    std::uniform_real_distribution<float> dist(0.0f, params_.global.noise);
    for (float &s : penalty) {
      s *= (1.0f + dist(rgen_));
    }
  }
  return penalty;
}

void GlobalPlacer::run() {
  runInitialLB();
  penalty_ = params_.global.penalty.initialValue;
  approximationDistance_ = initialApproximationDistance();
  penaltyCutoffDistance_ = initialPenaltyCutoffDistance();
  double nextPenaltyUpdateDistance = penaltyUpdateDistance();

  float lb = valueLB();
  float ub = std::numeric_limits<float>::infinity();
  int firstStep = params_.global.nbInitialSteps + 1;
  for (step_ = firstStep; step_ <= params_.global.maxNbSteps; ++step_) {
    std::cout << "#" << step_ << ":" << std::flush;
    runUB();
    ub = valueUB();
    std::cout << std::defaultfloat << std::setprecision(4) << "\tUB " << ub;

    float dist = leg_.meanDistance();
    std::cout << std::fixed << std::setprecision(1) << "\tDist " << dist / averageCellLength_;
    std::cout << std::flush;

    float gap = (ub - lb) / ub;
    // Stop if distance or the difference between LB and UB is small enough
    if (gap < params_.global.gapTolerance || dist < distanceTolerance()) {
      std::cout << std::endl;
      break;
    }
    if (dist < nextPenaltyUpdateDistance) {
      callback(PlacementStep::PenaltyUpdate, xPlacementUB_, yPlacementUB_);
      nextPenaltyUpdateDistance /= params_.global.penaltyUpdateBackoff;
    }
    for (int i = 0; i < params_.global.nbStepsBeforeRoughLegalization; ++i) {
      if (i != 0) {
        std::cout << "#" << step_ << ":\t........\t........";
        std::cout << std::flush;
      }
      runLB();
      lb = valueLB();
      std::cout << std::defaultfloat << std::setprecision(4) << "\tLB " << lb
                << std::endl;
    }
    penalty_ *= params_.global.penalty.updateFactor;
    penaltyCutoffDistance_ *= params_.global.penalty.cutoffDistanceUpdateFactor;
    approximationDistance_ *=
        params_.global.continuousModel.approximationDistanceUpdateFactor;
  }
  runUB();
  // De-densification only ever changes how cells were spread; the exported
  // placement (and everything downstream) must use the true cell sizes.
  restoreBaseWidths();
}

float GlobalPlacer::valueLB() const {
  return xtopo_.value(xPlacementLB_) + ytopo_.value(yPlacementLB_);
}

float GlobalPlacer::valueUB() const {
  return xtopo_.value(xPlacementUB_) + ytopo_.value(yPlacementUB_);
}

void GlobalPlacer::runInitialLB() {
  NetModel::Parameters params;
  params.netModel = params_.global.continuousModel.netModel;
  params.tolerance =
      params_.global.continuousModel.conjugateGradientErrorTolerance;
  params.maxNbIterations =
      params_.global.continuousModel.maxNbConjugateGradientSteps;
  g_dumpLabel = "step0-initLB-x";
  xPlacementLB_ = xtopo_.solveStar(params);
  g_dumpLabel = "step0-initLB-y";
  yPlacementLB_ = ytopo_.solveStar(params);
  std::cout << std::defaultfloat << std::setprecision(4) << "#0:\tLB "
            << valueLB() << std::endl;
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
  for (step_ = 1; step_ <= params_.global.nbInitialSteps; ++step_) {
    xPlacementLB_ = xtopo_.solve(xPlacementLB_, params);
    yPlacementLB_ = ytopo_.solve(yPlacementLB_, params);
    std::cout << std::defaultfloat << std::setprecision(4) << "#" << step_
              << ":\tLB " << valueLB() << std::endl;
    callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
  }
  // Simplify blending solutions by having a UB immediately
  xPlacementUB_ = xPlacementLB_;
  yPlacementUB_ = yPlacementLB_;
}

void GlobalPlacer::runLB() {
  // Compute the parameters for the continuous model solver
  NetModel::Parameters params;
  params.netModel = params_.global.continuousModel.netModel;
  params.approximationDistance = approximationDistance_;
  params.penaltyCutoffDistance = penaltyCutoffDistance_;
  params.tolerance =
      params_.global.continuousModel.conjugateGradientErrorTolerance;
  params.maxNbIterations =
      params_.global.continuousModel.maxNbConjugateGradientSteps;

  // Compute the per-cell penalty with randomization
  std::vector<float> penalty = computeIterationPerCellPenalty();

  float w = params_.global.penalty.targetBlending;
  std::vector<float> xTarget = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yTarget = blendPlacement(yPlacementLB_, yPlacementUB_, w);

  // Solve the continuous model (x and y independently).
  // Sequential (not async) so the dumped matrices are cleanly ordered.
  g_dumpLabel = "step" + std::to_string(step_) + "-LB-x";
  xPlacementLB_ = xtopo_.solveWithPenalty(xPlacementLB_, xTarget, penalty, params);
  g_dumpLabel = "step" + std::to_string(step_) + "-LB-y";
  yPlacementLB_ = ytopo_.solveWithPenalty(yPlacementLB_, yTarget, penalty, params);
  callback(PlacementStep::LowerBound, xPlacementLB_, yPlacementLB_);
}

void GlobalPlacer::runUB() {
  updateCellSizes();
  float w = params_.global.roughLegalization.targetBlending;
  std::vector<float> xTarget = blendPlacement(xPlacementLB_, xPlacementUB_, w);
  std::vector<float> yTarget = blendPlacement(yPlacementLB_, yPlacementUB_, w);
  leg_.updateCellTargetX(xTarget);
  leg_.updateCellTargetY(yTarget);
  leg_.run();
  xPlacementUB_ = leg_.spreadCoordX(xTarget);
  yPlacementUB_ = leg_.spreadCoordY(yTarget);
  callback(PlacementStep::UpperBound, xPlacementUB_, yPlacementUB_);
  applyDensification();
}

void GlobalPlacer::callback(PlacementStep step,
                            const std::vector<float> &xplace,
                            const std::vector<float> &yplace) {
  // --- instrumentation: append current coordinates to a CSV per step ---
  if (const char *dir = std::getenv("COLO_DUMP")) {
    static int coordSeq = 0;
    const char *name = step == PlacementStep::LowerBound ? "LB"
                       : step == PlacementStep::UpperBound ? "UB"
                       : step == PlacementStep::PenaltyUpdate ? "PEN"
                       : "DET";
    std::string path = std::string(dir) + "/coords.csv";
    bool header = (coordSeq == 0);
    std::ofstream os(path, std::ios::app);
    if (header) {
      os << "seq,step,type";
      for (int i = 0; i < (int)xplace.size(); ++i) os << ",x" << i << ",y" << i;
      os << "\n";
    }
    os << coordSeq++ << "," << step_ << "," << name;
    for (int i = 0; i < (int)xplace.size(); ++i)
      os << "," << xplace[i] << "," << yplace[i];
    os << "\n";
  }
  if (!callback_.has_value()) return;
  exportPlacement(circuit_, xplace, yplace);
  callback_.value()(step);
}

void GlobalPlacer::updateCellSizes() {
  if (circuit_.hasCellSizeUpdate_) {
    leg_.updateCellDemand(circuit_);
    averageCellLength_ = computeAverageCellSize();
    perCellPenalty_ = computePerCellPenalty();
    circuit_.hasCellSizeUpdate_ = false;
  }
}

void GlobalPlacer::updateNets() {
  if (circuit_.hasNetUpdate_) {
    xtopo_ = NetModel::xTopology(circuit_);
    ytopo_ = NetModel::yTopology(circuit_);
    circuit_.hasNetUpdate_ = false;
  }
}

void GlobalPlacer::applyDensification() {
  const DensificationParameters &dp = params_.global.densification;
  if (dp.mode == DensificationMode::Disabled) {
    return;
  }
  // Ramp the effect in gently over the first nbRampSteps iterations.
  float ramp = std::min(1.0f, (float)step_ / (float)dp.nbRampSteps);

  std::vector<float> factor = dp.mode == DensificationMode::Targeted
                                  ? targetedExpansion(ramp)
                                  : uniformExpansion(ramp);

  // Recompute widths from the baseline so factors never compound across steps.
  std::vector<int> widths = baseCellWidth_;
  for (int i = 0; i < circuit_.nbCells(); ++i) {
    if (circuit_.isFixed(i)) {
      continue;
    }
    int w = (int)std::lround(baseCellWidth_[i] * factor[i]);
    widths[i] = std::max(baseCellWidth_[i], w);
  }
  circuit_.setCellWidth(widths);  // trips hasCellSizeUpdate_ for the next UB
  densificationActive_ = true;
}

std::vector<float> GlobalPlacer::uniformExpansion(float ramp) const {
  const DensificationParameters &dp = params_.global.densification;
  // Largest factor that keeps total area under the target density, then cap.
  double feasible =
      baseDensity_ > 1e-9 ? dp.targetDensity / baseDensity_ : dp.maxFactor;
  double targetF = std::max(1.0, std::min(dp.maxFactor, feasible));
  float f = 1.0f + (float)(targetF - 1.0) * ramp;
  return std::vector<float>(circuit_.nbCells(), f);
}

std::vector<float> GlobalPlacer::targetedExpansion(float ramp) const {
  const DensificationParameters &dp = params_.global.densification;
  int n = circuit_.nbCells();
  std::vector<float> factor(n, 1.0f);
  if (baseDensity_ <= 1e-9) {
    return factor;
  }

  // Build a coarse density grid over the lower-bound (wirelength) placement,
  // where cells pile up: this is the congestion the router will fight.
  Rectangle area = leg_.placementArea();
  float w = (float)area.width();
  float h = (float)area.height();
  if (w <= 0.0f || h <= 0.0f) {
    return factor;
  }
  float binLen = std::max(1.0f, 4.0f * averageCellLength_);
  int nbx = std::min(64, std::max(1, (int)(w / binLen)));
  int nby = std::min(64, std::max(1, (int)(h / binLen)));
  float bw = w / nbx;
  float bh = h / nby;

  auto binIndex = [&](float pos, float lo, float sz, int nb) {
    int idx = (int)((pos - lo) / sz);
    return std::min(nb - 1, std::max(0, idx));
  };

  // Accumulate (true) cell area into the bin of each cell's LB position.
  std::vector<double> usage((size_t)nbx * nby, 0.0);
  std::vector<int> cellBin(n, -1);
  for (int i = 0; i < n; ++i) {
    if (circuit_.isFixed(i)) {
      continue;
    }
    int bx = binIndex(xPlacementLB_[i], area.minX, bw, nbx);
    int by = binIndex(yPlacementLB_[i], area.minY, bh, nby);
    int b = by * nbx + bx;
    cellBin[i] = b;
    usage[b] += (double)circuit_.area(i);
  }

  double binArea = (double)bw * (double)bh;
  for (int i = 0; i < n; ++i) {
    if (cellBin[i] < 0) {
      continue;
    }
    double localDensity = usage[cellBin[i]] / binArea;
    // How many times above the average density this region sits.
    double ratio = localDensity / baseDensity_;
    double over = std::max(0.0, ratio - 1.0);
    double f = 1.0 + dp.targetedStrength * over * ramp;
    factor[i] = (float)std::max(1.0, std::min((double)dp.maxFactor, f));
  }
  return factor;
}

void GlobalPlacer::restoreBaseWidths() {
  if (densificationActive_) {
    circuit_.setCellWidth(baseCellWidth_);
    densificationActive_ = false;
  }
}
}  // namespace coloquinte
