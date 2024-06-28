#include "ProcessedLAPPDFilter.h"

ProcessedLAPPDFilter::ProcessedLAPPDFilter() : Tool() {}

bool ProcessedLAPPDFilter::Initialise(std::string configfile, DataModel &data)
{

  /////////////////// Useful header ///////////////////////
  if (configfile != "")
    m_variables.Initialise(configfile); // loading config file
  // m_variables.Print();

  m_data = &data; // assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  FilterVerbosity = 0;
  m_variables.Get("FilterVerbosity", FilterVerbosity);
  requirePulsedAmp = 15;
  m_variables.Get("requirePulsedAmp", requirePulsedAmp);
  requirePulsedStripNumber = 7;
  m_variables.Get("requirePulsedStripNumber", requirePulsedStripNumber);

  saveAllLAPPDEvents = false;
  m_variables.Get("saveAllLAPPDEvents", saveAllLAPPDEvents);
  savePMTClusterEvents = false;
  m_variables.Get("savePMTClusterEvents", savePMTClusterEvents);

  MRDDataName = "FilteredMRDData";
  MRDNoVetoDataName = "FilteredMRDNoVetoData";
  m_variables.Get("MRDDataName", MRDDataName);
  m_variables.Get("MRDNoVetoDataName", MRDNoVetoDataName);

  AllLAPPDDataName = "FilteredAllLAPPDData";
  m_variables.Get("AllLAPPDDataName", AllLAPPDDataName);
  PMTClusterDataName = "FilteredPMTClusterData";
  m_variables.Get("PMTClusterDataName", PMTClusterDataName);

  filterType = "MRDtrack"; // MRDtrack and pmt cluster
  // filterType = "PMTCluster"; // pmt cluster only
  m_variables.Get("filterType", filterType);

  FilteredMRD = new BoostStore(false, 2);        // all events with MRD in it //if gotEventMRD = true;
  FilteredMRDNoVeto = new BoostStore(false, 2);  // all events with MRD but no veto// if gotEventMRD && gotEventMRDNoVeto = true
  FilteredAllLAPPD = new BoostStore(false, 2);   // all events with LAPPD data
  FilteredPMTCluster = new BoostStore(false, 2); // all events with PMT cluster

  gotEventMRD = false;
  gotEventMRDNoVeto = false;
  gotEventPMTCluster = false;

  EventPMTClusterNumber = 0;
  EventMRDNumber = 0;
  EventMRDNoVetoNumber = 0;
  CurrentExeNumber = 0;
  pairedEventNumber = 0;
  PsecDataNumber = 0;
  PMTClusterEventNum = 0;

  return true;
}

bool ProcessedLAPPDFilter::Execute()
{

  m_data->Stores.at("ANNIEEvent")->Get("DataStreams", DataStreams);

  if (FilterVerbosity > 2)
  {
    for (auto it = DataStreams.begin(); it != DataStreams.end(); ++it)
    {
      cout << "DataStream: " << it->first << " " << it->second << ", ";
    }
    cout << endl;
  }

  if (filterType == "MRDtrack")
  {
    gotEventMRD = false;
    gotEventPMTCluster = false;

    std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;
    bool get_clusters = m_data->CStore.Get("ClusterMap", m_all_clusters);
    if (!get_clusters)
    {
      std::cout << "ProcessedLAPPDFilter tool: No clusters found!" << std::endl;
      return false;
    }
    if ((int)m_all_clusters->size() > 0)
      gotEventPMTCluster = true;

    if (savePMTClusterEvents && gotEventPMTCluster)
    {
      GotANNIEEventAndSave(FilteredPMTCluster, PMTClusterDataName);
      PMTClusterEventNum += 1;
    }

    if (DataStreams["LAPPD"] == false)
    {
      return true;
    }

    CurrentExeNumber += 1;

    std::map<uint64_t, uint64_t> LAPPDBeamgate_ns;
    m_data->Stores["ANNIEEvent"]->Get("LAPPDBeamgate_ns", LAPPDBeamgate_ns);
    PsecDataNumber += LAPPDBeamgate_ns.size();
    m_data->Stores["ANNIEEvent"]->Set("LAPPDBeamgate_ns", LAPPDBeamgate_ns);

    std::vector<std::vector<int>> MrdTimeClusters;
    bool get_clusters_mrd = m_data->CStore.Get("MrdTimeClusters", MrdTimeClusters);
    if (!get_clusters_mrd)
    {
      std::cout << "ProcessedLAPPDFilter tool: No MRD clusters found! Did you run the TimeClustering tool?" << std::endl;
      return false;
    }
    if (MrdTimeClusters.size() > 0)
      gotEventMRD = true;

    bool gotStoreNoVeto = m_data->Stores.at("RecoEvent")->Get("NoVeto", gotEventMRDNoVeto);
    if (!gotStoreNoVeto)
    {
      cout << "The Stage1 Filter doesn't work because it need the Veto information form EventSelector Tool. Please include that." << endl;
      return false;
    }

    if (saveAllLAPPDEvents)
    {
      GotANNIEEventAndSave(FilteredAllLAPPD, AllLAPPDDataName);
    }

    if (gotEventPMTCluster)
    {
      EventPMTClusterNumber += 1;
      if (FilterVerbosity > 0)
        cout << "Got an event with PMT cluster, got event number " << EventPMTClusterNumber << endl;
      if (gotEventMRD)
      {
        GotANNIEEventAndSave(FilteredMRD, MRDDataName);
        EventMRDNumber += 1;
        if (FilterVerbosity > 0)
          cout << "Got an event with PMT cluster, MRD hits, got event number " << EventMRDNumber << endl;
        if (gotEventMRDNoVeto)
        {
          GotANNIEEventAndSave(FilteredMRDNoVeto, MRDNoVetoDataName);
          EventMRDNoVetoNumber += 1;
          if (FilterVerbosity > 0)
            cout << "Got an event with PMT cluster, MRD hits and no veto, got event number " << EventMRDNoVetoNumber << endl;
        }
      }

      if (CurrentExeNumber % 10 == 0)
      {
        cout << "Filter event number: " << CurrentExeNumber << ", PsecDataNumber " << PsecDataNumber << ", LAPPD events with PMT clusters: " << EventPMTClusterNumber << ", also with MRD tracks: " << EventMRDNumber << ", also with MRD tracks and no veto: " << EventMRDNoVetoNumber << ", events with PMT cluster (may not have LAPPD): " << PMTClusterEventNum << endl;
      }

      m_data->Stores["ANNIEEvent"]->Get("RunNumber", currentRunNumber);
    }
  }

  if (filterType == "PMTCluster")
  {
    CurrentExeNumber += 1;
    std::map<uint64_t, uint64_t> LAPPDBeamgate_ns;
    m_data->Stores["ANNIEEvent"]->Get("LAPPDBeamgate_ns", LAPPDBeamgate_ns);
    PsecDataNumber += LAPPDBeamgate_ns.size();

    gotEventPMTCluster = false;

    std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;
    bool get_clusters = m_data->CStore.Get("ClusterMap", m_all_clusters);
    if (!get_clusters)
    {
      std::cout << "ProcessedLAPPDFilter tool: No clusters found!" << std::endl;
      return false;
    }
    if ((int)m_all_clusters->size() > 0)
      gotEventPMTCluster = true;

    if (gotEventPMTCluster)
    {
      EventPMTClusterNumber += 1;
      if (FilterVerbosity > 0)
        cout << "Got an event with PMT cluster, got event number " << EventPMTClusterNumber << endl;
      GotANNIEEventAndSave(FilteredMRD, MRDDataName);
      if (CurrentExeNumber % 50 == 0)
      {
        cout << "Filter event number: " << CurrentExeNumber << ", PsecDataNumber " << PsecDataNumber << ", events with PMT clusters: " << EventPMTClusterNumber << endl;
      }

      m_data->Stores["ANNIEEvent"]->Get("RunNumber", currentRunNumber);
    }
  }

  /*
    if (filterType == "LAPPDPulsedStripNumber")
    {
      CurrentExeNumber += 1;
      int pulsedStripNumber = 0;
      std::vector<std::vector<double>> LAPPDVectorSide0 = std::vector<std::vector<double>>(30, std::vector<double>(256));
      // loop the vector. For each strip, loop and find maximum, if the maximum > 15, then pulsedStripNumber += 1
      m_data->Stores["ANNIEEvent"]->Get("LAPPDSigVecSide0", LAPPDVectorSide0);
      for (int i = 0; i < 30; i++)
      {
        double max = 0;
        for (int j = 0; j < 256; j++)
        {
          if (LAPPDVectorSide0[i][j] > max)
            max = LAPPDVectorSide0[i][j];
          if (FilterVerbosity > requirePulsedStripNumber)
            cout << "LAPPDVectorSide0[" << i << "][" << j << "] = " << LAPPDVectorSide0[i][j] << endl;
        }
        if (max > requirePulsedAmp)
          pulsedStripNumber += 1;
      }
      if (pulsedStripNumber > requirePulsedStripNumber)
      {
        GotANNIEEventAndSave(FilteredMRD, MRDDataName);
        EventMRDNumber += 1;
      }

      if (CurrentExeNumber % 50 == 0)
      {
        cout << "Filter event number " << CurrentExeNumber << ", pulsedStripNumber more than " << requirePulsedStripNumber << " events " << EventMRDNumber << endl;
      }
    }
  */

  return true;
}

bool ProcessedLAPPDFilter::Finalise()
{
  FilteredMRD->Close();
  FilteredMRDNoVeto->Close();
  FilteredAllLAPPD->Close();
  FilteredPMTCluster->Close();

  FilteredMRD->Delete();
  FilteredMRDNoVeto->Delete();
  FilteredAllLAPPD->Delete();
  FilteredPMTCluster->Delete();

  if (filterType == "MRDtrack")
  {
    cout << "Current Run " << currentRunNumber << endl;
    cout << "Filter got " << CurrentExeNumber << " events in total" << endl; // this is the total number of events
    cout << "Got " << PsecDataNumber << " psec data objects on all LAPPD" << endl;
    cout << "Got " << EventPMTClusterNumber << " events with PMT clusters" << endl;
    cout << "Got " << EventMRDNumber << " events with PMT clusters and MRD hits" << endl;
    cout << "Got " << EventMRDNoVetoNumber << " events with PMT clusters and MRD hits and no veto" << endl;
    cout << "Got " << PMTClusterEventNum << " events with PMT clusters and saved" << endl;
    cout << "Saved " << MRDDataName << " successfully" << endl;
    cout << "Saved " << MRDNoVetoDataName << " successfully" << endl;
    cout << "Saved " << AllLAPPDDataName << " successfully" << endl;
    cout << "Saved " << PMTClusterDataName << " successfully" << endl;

    std::ofstream file("FilterStat.txt", std::ios::app);
    file << currentRunNumber << " " << CurrentExeNumber << " " << EventPMTClusterNumber << " " << EventMRDNumber << " " << EventMRDNoVetoNumber << " "
         << "\n";
    file.close();
  }
  if (filterType == "PMTCluster")
  {
    cout << "Current Run " << currentRunNumber << endl;
    cout << "Filter got " << CurrentExeNumber << " events in total" << endl; // this is the total number of events
    cout << "Got " << EventPMTClusterNumber << " events with PMT clusters" << endl;
    cout << "Saved " << MRDDataName << " successfully" << endl;
  }
  if (filterType == "LAPPDTiming")
  {
    cout << "Filter got " << CurrentExeNumber << " events in total" << endl; // this is the total number of events
    cout << "Got " << EventMRDNumber << " events in beamgate window" << endl;
    cout << "Saved " << MRDDataName << " successfully" << endl;
  }

  if (filterType == "LAPPDPulsedStripNumber")
  {
    cout << "Filter got " << CurrentExeNumber << " events in total" << endl; // this is the total number of events
    cout << "Got " << EventMRDNumber << " events with more than 5 pulsed strips" << endl;
    cout << "Saved " << MRDDataName << " successfully" << endl;
  }

  return true;
}

bool ProcessedLAPPDFilter::GotANNIEEventAndSave(BoostStore *BS, string savePath)
{
  std::map<unsigned long, std::vector<Hit>> *AuxHitsOriginal = nullptr;
  std::map<unsigned long, std::vector<Hit>> *HitsOriginal = nullptr;
  m_data->Stores["ANNIEEvent"]->Get("AuxHits", AuxHitsOriginal);
  m_data->Stores["ANNIEEvent"]->Get("Hits", HitsOriginal);
  auto AuxHitsCopy = new std::map<unsigned long, std::vector<Hit>>(*AuxHitsOriginal);
  auto HitsCopy = new std::map<unsigned long, std::vector<Hit>>(*HitsOriginal);
  BS->Set("AuxHits", AuxHitsCopy);
  BS->Set("Hits", HitsCopy);

  std::map<unsigned long, std::vector<Hit>> TDCData;
  BeamStatus BeamStatus;
  uint64_t CTCTimestamp, RunStartTime;
  std::map<std::string, bool> DataStreams;
  uint32_t EventNumber, TriggerWord;
  uint64_t EventTimeTank;
  TimeClass EventTimeMRD;
  std::map<std::string, int> MRDLoopbackTDC;
  std::string MRDTriggerType;
  std::map<unsigned long, std::vector<int>> RawAcqSize;
  std::map<unsigned long, std::vector<std::vector<ADCPulse>>> RecoADCData, RecoAuxADCData;
  int PartNumber, RunNumber, RunType, SubrunNumber, TriggerExtended;
  TriggerClass TriggerData;
  m_data->Stores["ANNIEEvent"]->Get("TDCData", TDCData);
  BS->Set("TDCData", TDCData);
  m_data->Stores["ANNIEEvent"]->Get("DataStreams", DataStreams);
  BS->Set("DataStreams", DataStreams);
  m_data->Stores["ANNIEEvent"]->Get("BeamStatus", BeamStatus);
  BS->Set("BeamStatus", BeamStatus);
  m_data->Stores["ANNIEEvent"]->Get("RunStartTime", RunStartTime);
  BS->Set("RunStartTime", RunStartTime);
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", EventNumber);
  BS->Set("EventNumber", EventNumber);
  m_data->Stores["ANNIEEvent"]->Get("TriggerWord", TriggerWord);
  BS->Set("TriggerWord", TriggerWord);
  m_data->Stores["ANNIEEvent"]->Get("EventTimeMRD", EventTimeMRD);
  BS->Set("EventTimeMRD", EventTimeMRD);
  m_data->Stores["ANNIEEvent"]->Get("EventTimeTank", EventTimeTank);
  BS->Set("EventTimeTank", EventTimeTank);
  m_data->Stores["ANNIEEvent"]->Get("MRDLoopbackTDC", MRDLoopbackTDC);
  BS->Set("MRDLoopbackTDC", MRDLoopbackTDC);
  m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType", MRDTriggerType);
  BS->Set("MRDTriggerType", MRDTriggerType);
  m_data->Stores["ANNIEEvent"]->Get("RawAcqSize", RawAcqSize);
  BS->Set("RawAcqSize", RawAcqSize);
  m_data->Stores["ANNIEEvent"]->Get("RecoADCData", RecoADCData);
  BS->Set("RecoADCData", RecoADCData);
  m_data->Stores["ANNIEEvent"]->Get("RecoAuxADCData", RecoAuxADCData);
  BS->Set("RecoAuxADCData", RecoAuxADCData);
  m_data->Stores["ANNIEEvent"]->Get("PartNumber", PartNumber);
  BS->Set("PartNumber", PartNumber);
  m_data->Stores["ANNIEEvent"]->Get("RunNumber", RunNumber);
  BS->Set("RunNumber", RunNumber);
  m_data->Stores["ANNIEEvent"]->Get("RunType", RunType);
  BS->Set("RunType", RunType);
  m_data->Stores["ANNIEEvent"]->Get("SubrunNumber", SubrunNumber);
  BS->Set("SubrunNumber", SubrunNumber);
  m_data->Stores["ANNIEEvent"]->Get("TriggerExtended", TriggerExtended);
  BS->Set("TriggerExtended", TriggerExtended);
  m_data->Stores["ANNIEEvent"]->Get("TriggerData", TriggerData);
  BS->Set("TriggerData", TriggerData);
  m_data->Stores["ANNIEEvent"]->Get("CTCTimestamp", CTCTimestamp);
  BS->Set("CTCTimestamp", CTCTimestamp);

  std::map<uint64_t, uint32_t> GroupedTrigger;
  m_data->Stores["ANNIEEvent"]->Get("GroupedTrigger", GroupedTrigger);
  BS->Set("GroupedTrigger", GroupedTrigger);

  int PrimaryTriggerWord;
  m_data->Stores["ANNIEEvent"]->Get("PrimaryTriggerWord", PrimaryTriggerWord);
  BS->Set("PrimaryTriggerWord", PrimaryTriggerWord);

  uint64_t primaryTrigTime;
  m_data->Stores["ANNIEEvent"]->Get("PrimaryTriggerTime", primaryTrigTime);
  BS->Set("PrimaryTriggerTime", primaryTrigTime);

  bool NCExtended, CCExtended;
  m_data->Stores["ANNIEEvent"]->Get("NCExtended", NCExtended);
  BS->Set("NCExtended", NCExtended);
  m_data->Stores["ANNIEEvent"]->Get("CCExtended", CCExtended);
  BS->Set("CCExtended", CCExtended);

  string TriggerType;
  m_data->Stores["ANNIEEvent"]->Get("TriggerType", TriggerType);
  BS->Set("TriggerType", TriggerType);

  int CTCWordExtended;
  m_data->Stores["ANNIEEvent"]->Get("CTCWordExtended", CTCWordExtended);
  BS->Set("CTCWordExtended", CTCWordExtended);

  std::map<uint64_t, PsecData> LAPPDDataMap;
  std::map<uint64_t, uint64_t> LAPPDBeamgate_ns;
  std::map<uint64_t, uint64_t> LAPPDTimeStamps_ns; // data and key are the same
  std::map<uint64_t, uint64_t> LAPPDTimeStampsRaw;
  std::map<uint64_t, uint64_t> LAPPDBeamgatesRaw;
  std::map<uint64_t, uint64_t> LAPPDOffsets;
  std::map<uint64_t, int> LAPPDTSCorrection;
  std::map<uint64_t, int> LAPPDBGCorrection;
  std::map<uint64_t, int> LAPPDOSInMinusPS;

  m_data->Stores["ANNIEEvent"]->Get("LAPPDDataMap", LAPPDDataMap);
  BS->Set("LAPPDDataMap", LAPPDDataMap);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDBeamgate_ns", LAPPDBeamgate_ns);
  BS->Set("LAPPDBeamgate_ns", LAPPDBeamgate_ns);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDTimeStamps_ns", LAPPDTimeStamps_ns);
  BS->Set("LAPPDTimeStamps_ns", LAPPDTimeStamps_ns);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDTimeStampsRaw", LAPPDTimeStampsRaw);
  BS->Set("LAPPDTimeStampsRaw", LAPPDTimeStampsRaw);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDBeamgatesRaw", LAPPDBeamgatesRaw);
  BS->Set("LAPPDBeamgatesRaw", LAPPDBeamgatesRaw);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDOffsets", LAPPDOffsets);
  BS->Set("LAPPDOffsets", LAPPDOffsets);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDTSCorrection", LAPPDTSCorrection);
  BS->Set("LAPPDTSCorrection", LAPPDTSCorrection);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDBGCorrection", LAPPDBGCorrection);
  BS->Set("LAPPDBGCorrection", LAPPDBGCorrection);
  m_data->Stores["ANNIEEvent"]->Get("LAPPDOSInMinusPS", LAPPDOSInMinusPS);
  BS->Set("LAPPDOSInMinusPS", LAPPDOSInMinusPS);

  uint64_t beamInfoTime = 0;
  int64_t timeDiff = -9999;
  double defaultVal = -9999.;
  int beam_good = 0;
  double E_TOR860, E_TOR875, THCURR, BTJT2, HP875, VP875, HPTG1, VPTG1, HPTG2, VPTG2, BTH2T2;
  m_data->Stores["ANNIEEvent"]->Get("BeamInfoTime", beamInfoTime);
  BS->Set("BeamInfoTime", beamInfoTime);
  m_data->Stores["ANNIEEvent"]->Get("BeamInfoTimeToTriggerDiff", timeDiff);
  BS->Set("BeamInfoTimeToTriggerDiff", timeDiff);
  m_data->Stores["ANNIEEvent"]->Get("beam_E_TOR860", E_TOR860);
  BS->Set("beam_E_TOR860", E_TOR860);
  m_data->Stores["ANNIEEvent"]->Get("beam_E_TOR875", E_TOR875);
  BS->Set("beam_E_TOR875", E_TOR875);
  m_data->Stores["ANNIEEvent"]->Get("beam_THCURR", THCURR);
  BS->Set("beam_THCURR", THCURR);
  m_data->Stores["ANNIEEvent"]->Get("beam_BTJT2", BTJT2);
  BS->Set("beam_BTJT2", BTJT2);
  m_data->Stores["ANNIEEvent"]->Get("beam_HP875", HP875);
  BS->Set("beam_HP875", HP875);
  m_data->Stores["ANNIEEvent"]->Get("beam_VP875", VP875);
  BS->Set("beam_VP875", VP875);
  m_data->Stores["ANNIEEvent"]->Get("beam_HPTG1", HPTG1);
  BS->Set("beam_HPTG1", HPTG1);
  m_data->Stores["ANNIEEvent"]->Get("beam_VPTG1", VPTG1);
  BS->Set("beam_VPTG1", VPTG1);
  m_data->Stores["ANNIEEvent"]->Get("beam_HPTG2", HPTG2);
  BS->Set("beam_HPTG2", HPTG2);
  m_data->Stores["ANNIEEvent"]->Get("beam_VPTG2", VPTG2);
  BS->Set("beam_VPTG2", VPTG2);
  m_data->Stores["ANNIEEvent"]->Get("beam_BTH2T2", BTH2T2);
  BS->Set("beam_BTH2T2", BTH2T2);
  m_data->Stores["ANNIEEvent"]->Get("beam_good", beam_good);
  BS->Set("beam_good", beam_good);

  BS->Save(savePath);
  if (FilterVerbosity > 2)
    cout << "Saved to " << savePath << " successfully" << endl;
  BS->Delete();

  // removd and clean the data pointers from here
  /* std::map<unsigned long, std::vector<Hit>> *AuxHitsOriginal = nullptr;
 std::map<unsigned long, std::vector<Hit>> *HitsOriginal = nullptr;
 m_data->Stores["ANNIEEvent"]->Get("AuxHits", AuxHitsOriginal);
 m_data->Stores["ANNIEEvent"]->Get("Hits", HitsOriginal);
 auto AuxHitsCopy = new std::map<unsigned long, std::vector<Hit>>(*AuxHitsOriginal);
 auto HitsCopy = new std::map<unsigned long, std::vector<Hit>>(*HitsOriginal);
 BS->Set("AuxHits", AuxHitsCopy);
 BS->Set("Hits", HitsCopy);

  AuxHitsCopy->clear();
  HitsCopy->clear();
  delete AuxHitsCopy;
  delete HitsCopy;*/

  return true;
}