void inspect_scan_file(const char* fname) {
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) { printf("BADFILE\n"); return; }
    TList* keys = f->GetListOfKeys();
    printf("NKEYS|%d\n", keys->GetEntries());
    for (int i = 0; i < keys->GetEntries(); i++) {
        TKey* k = (TKey*)keys->At(i);
        TString name = k->GetName();
        if (name == "h_eventsBeforeSelection" || name == "h_eventsAfterSelection") {
            TObject* obj = k->ReadObj();
            if (obj->InheritsFrom("TH1")) {
                TH1* h = (TH1*)obj;
                printf("HIST|%s|%f\n", name.Data(), h->GetEntries());
            }
        }
    }
    f->Close();
}
