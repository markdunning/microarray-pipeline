produceGEOSampleInfoTemplate <-
function(lumiNormalized, lib.mapping=NULL, fileName='GEOsampleInfo.txt') {	
	
	link <- "http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html"
	templateTitle <- c("Sample_title", "Sample_channel_count","Sample_source_name_ch1","Sample_organism_ch1", "Sample_characteristics_ch1", "Sample_molecule_ch1","Sample_extract_protocol_ch1","Sample_label_ch1", "Sample_label_protocol_ch1", "Sample_hyb_protocol","Sample_scan_protocol","Sample_description","Sample_data_processing","Sample_platform_id", "Sample_supplementary_file")
	if (is(lumiNormalized, 'ExpressionSetIllumina')) {

		labels <- as.character(pData(lumiNormalized)[,"Sample_Name"])

#		chipInfo <- getChipInfo(lumiNormalized, lib.mapping=lib.mapping)

#		chipVersion <- chipInfo$chipVersion[1]

#		chipVersion <- sub("(_V[0-9.]*)_.*$", "\\1", chipVersion)

		chipVersion <- annotation(lumiNormalized)

		organism <- strsplit(annotation(lumiNormalized), "v")[[1]][1]

		organism <- switch(organism,
			'Rat'='Rattus norvegicus',
			'Human'="Homo sapiens",
			'Mouse'='Mus musculus')


		templateContent <- c("","1","",organism,"","total RNA","standard as recommended by illumina","Cy3","standard as recommended by illumina","standard as recommended by illumina","standard as recommended by illumina","","",chipVersion, "none")

	} else if (is(lumiNormalized, 'MethyLumiM')) {
		labels <- as.character(sampleNames(lumiNormalized))
		chipVersion <- 'unknown'
		organism <- 'unknown'
		templateContent <- c("","1","",organism,"","genomic DNA","standard as recommended by illumina","Cy3","standard as recommended by illumina","standard as recommended by illumina","standard as recommended by illumina","","",chipVersion,"none")
	} else if (is(lumiNormalized, 'matrix') || is(lumiNormalized, 'ExpressionSet')) {
		if (is(lumiNormalized, 'ExpressionSet'))  lumiNormalized <- exprs(lumiNormalized)
		labels <- colnames(lumiNormalized)
		chipVersion <- 'unknown'
		organism <- 'unknown'
		templateContent <- c("","1","",organism,"","","","Cy3","","standard as recommended by manufacturer","standard as recommended by manufacturer","","", "","none")
	} else {
		cat("The input object should be an object of LumiBatch, MethyLumiM, matrix or other ExpressionSet inherited class!\n")
	}
	
	## add code of parsing processing history 

	 preprocessMethod <- 'The bead-level data were preprocessed using BASH (Cairns et al. (2008) Bioinformatics 24(24):2921-2), a function from the beadarray package (Dunning et al (2007) Bioinformatics 23(16):2183-4) in Bioconductor, and also log base 2 transformed and quantile normalised using the beadarray packagee'		

	templateContent[templateTitle == "Sample_data_processing"] <- preprocessMethod
	templateContent[templateTitle == "Sample_characteristics_ch1"] <- as.character(pData(lumiNormalized)[,"Sample_Group"])

	template <- templateTitle
	for (i in seq(labels)) {
		template <- rbind(template, templateContent)
	}


	template <- cbind(c('sampleID', labels), template)
	template <- cbind(c('Sample_title', as.character(pData(lumiNormalized)[,"Sample_Name"])), template)

	template
	colnames(template) <- template[1,]
	template <- template[-1,]

	if (!is.null(fileName)) {
		cat('# For the detailed definition of the column names, please refer to ', link, '\n', sep='', file=fileName)
		write.table(template, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)		
	} else {
		return(template)
	}

      

}
