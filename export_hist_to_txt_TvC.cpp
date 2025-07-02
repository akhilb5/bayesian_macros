void export_hist_to_txt_TvC(TH2* hist, const char* output_filename = "histogram_output.txt") {
    if (!hist) {
        std::cerr << "Histogram is null." << std::endl;
        return;
    }

    std::ofstream txt(output_filename);
    if (!txt.is_open()) {
        std::cerr << "Could not open file: " << output_filename << std::endl;
        return;
    }

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();
    double BinContent = 0.0;
    int count = 0;
    for (int xi = 1; xi <= nBinsX; ++xi) {
        for(int yi = 1; yi <= nBinsY; ++yi) {
            BinContent = hist->GetBinContent(xi, yi);
            //if(BinContent>0.0)
            txt << xi << "\t" << yi << "\t" << BinContent << "\n";
            count++;
        }
    }
    std::cout << "Wrote " << count << std::endl;
    // for (int bin = 1; bin <= ; ++bin) {
    //     double center = hist->GetBinCenter(bin);
    //     double content = hist->GetBinContent(bin);
    //     //double error = hist->GetBinError(bin);
    //     txt << center << "\t" << content << "\n";
    // }

    txt.close();
    std::cout << "Histogram exported to " << output_filename << std::endl;
}
