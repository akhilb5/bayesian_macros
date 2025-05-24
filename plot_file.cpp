void plotHistogram() {
    // Open the file
    TTree *tree = new TTree("tree", "data from txt");
    tree->ReadFile("Level1a.txt", "x/D:y/D");

    // Create a histogram for y values
    TH1D *hist = new TH1D("hist", "Histogram of Y;Y Value;Frequency", 50, 0, 7000);

    // Fill histogram with 'y' values from the tree
    tree->Draw("y >> hist");

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    hist->Draw();
}
s[0] = 0.0873984, s[1] = 0.402048, s[2] = 0.51055


