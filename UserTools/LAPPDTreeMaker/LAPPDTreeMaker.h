#ifndef LAPPDTreeMaker_H
#define LAPPDTreeMaker_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "LAPPDPulse.h"
#include "LAPPDHit.h"
#include "Geometry.h"
#include "Position.h"

/**
 * \class LAPPDTreeMaker
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
 *
 * $Author: Yue Feng $
 * $Date: 2024/02/03 16:10:00 $
 * Contact: yuef@iastate.edu
 */
class LAPPDTreeMaker : public Tool
{

public:
    LAPPDTreeMaker();                                         ///< Simple constructor
    bool Initialise(std::string configfile, DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
    bool Execute();                                           ///< Execute function used to perform Tool purpose.
    bool Finalise();                                          ///< Finalise function used to clean up resources.
    void CleanVariables();
    bool FillPulseTree();
    bool FillHitTree();
    bool FillWaveformTree();
    bool FillPPSTimestamp();
    bool FillLAPPDDataTimeStamp();
    bool LoadRunInfoFromRaw();
    bool LoadRunInfoFromANNIEEvent();
    bool GroupTriggerByBeam();
    bool GroupTriggerByLaser();
    bool GroupPPSTrigger();
    bool FillTriggerTree();
    bool FillGroupedTriggerTree();
    void CleanTriggers();

private:
    TFile *file;
    TTree *fPulse;
    TTree *fHit;
    TTree *fWaveform;
    TTree *fTimeStamp;
    TTree *fTrigger;
    TTree *fGroupedTrigger;

    int treeMakerVerbosity;
    bool LAPPDana;
    bool LoadingPPS;

    bool LoadPulse;
    bool LoadHit;
    bool LoadWaveform;
    bool LoadLAPPDDataTimeStamp;
    bool LoadPPSTimestamp;
    bool LoadRunInfoRaw;
    bool LoadRunInfoANNIEEvent;
    bool LoadTriggerInfo;
    bool LoadGroupedTriggerInfo;
    string LoadGroupOption; // using beam, or laser, or other.

    std::map<unsigned long, vector<Waveform<double>>> lappdData;
    std::map<unsigned long, vector<double>> waveformMax; // strip number+30*side, value
    std::map<unsigned long, vector<double>> waveformRMS;
    std::map<unsigned long, vector<double>> waveformMaxLast;
    std::map<unsigned long, vector<double>> waveformMaxNearing;
    std::map<unsigned long, vector<int>> waveformMaxTimeBin;

    map<uint64_t, std::vector<uint32_t>> *TriggerWordMap;
    unsigned long CTCTimeStamp;
    std::vector<uint32_t> CTCTriggerWord;
    int trigNumInThisMap;
    int trigIndexInThisMap;
    int TriggerGroupNumInThisEvent;
    int groupedTriggerType;

    vector<uint32_t> unGroupedTriggerWords;
    vector<unsigned long> unGroupedTriggerTimestamps;
    vector<vector<uint32_t>> groupedTriggerWordsVector;
    vector<vector<unsigned long>> groupedTriggerTimestampsVector;
    vector<int> groupedTriggerByType;
    vector<uint32_t> groupedTriggerWords;
    vector<unsigned long> groupedTriggerTimestamps;

    double waveformMaxValue;
    double waveformRMSValue;
    bool waveformMaxFoundNear;
    double waveformMaxNearingValue;
    int waveformMaxTimeBinValue;

    std::map<unsigned long, vector<vector<LAPPDPulse>>> lappdPulses;
    std::map<unsigned long, vector<LAPPDHit>> lappdHits;
    Geometry *_geom;
    string treeMakerInputPulseLabel;
    string treeMakerInputHitLabel;
    string treeMakerOutputFileName;

    int EventNumber; // this is for LAPPD relate event number, data or PPS, not trigger!!
    int StripNumber;

    int LAPPD_ID;
    int ChannelID;
    double PeakTime;
    double Charge;
    double PeakAmp;
    double PulseStart;
    double PulseEnd;
    double PulseSize;
    int PulseSide;
    double PulseThreshold;
    double PulseBaseline;

    double HitTime;
    double HitAmp;
    double XPosTank;
    double YPosTank;
    double ZPosTank;
    double ParallelPos;
    double TransversePos;
    double Pulse1StartTime;
    double Pulse2StartTime;
    double Pulse1LastTime;
    double Pulse2LastTime;

    unsigned long LAPPDDataTimeStampUL;
    unsigned long LAPPDDataBeamgateUL;
    unsigned long LAPPDDataTimestampPart1;
    unsigned long LAPPDDataBeamgatePart1;
    double LAPPDDataTimestampPart2;
    double LAPPDDataBeamgatePart2;

    vector<unsigned long> pps_vector;
    vector<unsigned long> pps_count_vector;
    Long64_t ppsDiff;
    unsigned long ppsCount0;
    unsigned long ppsCount1;
    unsigned long ppsTime0; // in clock tick, without 3.125ns
    unsigned long ppsTime1; // in clock tick, without 3.125ns

    int RunNumber;
    int SubRunNumber;
    int PartFileNumber;
};

#endif