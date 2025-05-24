void export_hist_to_txt(TH1* hist, const char* output_filename = "histogram_output.txt") {
    if (!hist) {
        std::cerr << "Histogram is null." << std::endl;
        return;
    }

    std::ofstream txt(output_filename);
    if (!txt.is_open()) {
        std::cerr << "Could not open file: " << output_filename << std::endl;
        return;
    }

    // Write histogram bin center and content (optionally also include errors)
    //txt << "# BinCenter\tContent\tError\n";
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        double center = hist->GetBinCenter(bin);
        double content = hist->GetBinContent(bin);
        //double error = hist->GetBinError(bin);
        txt << center << "\t" << content << "\n";
    }

    txt.close();
    std::cout << "Histogram exported to " << output_filename << std::endl;
}
