#include "Analysis/DrellYan/interface/Sample.h"
#include "Analysis/DrellYan/interface/Yield.h"
#ifdef __MAKECINT__ 

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace dy;
#pragma link C++ typedef dy::YieldMap;
#pragma link C++ typedef dy::SampleMap;
#pragma link C++ enum dy::Sample::value_type;
#pragma link C++ class dy::Sample+;
#pragma link C++ class dy::Sample::Info+;
#pragma link C++ class dy::Yield+;
#pragma link C++ class dy::Yield::value_t+;
#pragma link C++ class std::map<dy::Sample::value_type, dy::Yield>+;
#pragma link C++ class std::map<dy::Sample::value_type, dy::Sample::Info>+;
#pragma link C++ function dy::GetSampleFromName;
#pragma link C++ function dy::GetSampleFromNumber;
#pragma link C++ function dy::IsSample;
#pragma link C++ function dy::GetSampleInfo;
#pragma link C++ function dy::GetSampleInfo;
#pragma link C++ function dy::GetSampleInfo;
#pragma link C++ function dy::GetSampleTChain;
#pragma link C++ function dy::GetSampleMap;
#pragma link C++ function dy::GetYieldFromHist;
#pragma link C++ function dy::GetRecoYieldFromLabel;
#pragma link C++ function dy::GetGenYieldFromLabel;
#pragma link C++ function dy::GetYieldString;
#pragma link C++ function dy::GetRecoYieldMap;
#pragma link C++ function dy::GetGenYieldMap;

#endif

