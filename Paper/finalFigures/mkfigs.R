library(pdftools)

# combine PDFs, the file names will already be in order
pdf_combine(list.files("Paper/finalFigures/Rasterized", pattern = "(F|S)[0-9]", full.names = TRUE), output = "Paper/finalFigures/Figures_raster.pdf")


#pdf_combine(list.files("Paper/finalFigures/Rasterized/With_legends", pattern = "S[0-9]", full.names = TRUE), output = "Paper/finalFigures/Figures_supplemental_w_legends.pdf")
