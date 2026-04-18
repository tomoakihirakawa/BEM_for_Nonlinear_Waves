#pragma once

#include "Network.hpp"

void dumpDebugCornerPointState(Network* water, const char* stage, int time_step, int rk_step);
void installCrashBacktraceIfRequested();
void enableRealFieldFmmDefaultForTimeDomain();
void logCornerConnectedNeumannLinesAfterBoundaryTypes(const Network* water, const char* phase);
