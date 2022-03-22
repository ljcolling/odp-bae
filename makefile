all: analysis, document
	make analysis
	make document

analysis : 
	First get the data
	@echo "Downloading the data... this may take a while"
	@curl https://files.osf.io/v1/resources/dbgh6/providers/osfstorage/\?zip\= -o data.zip
	@echo "Unizipping the data"
	@unzip ./data.zip -d data
	@echo "Moving code files into data folder"
	@cp ./code/matlab ./data/
	@cp ./code/do_analysis.m ./data
	@cp ./code/SVM_ECOC_ERP_Decoding_PreviousTrial_1.m './data/Experiment 1 data and analysis code/'
	@cp ./code/SVM_ECOC_ERP_Decoding_PreviousTrial_2.m './data/Experiment 2 data and analysis code/'
	@echo "Now running the analysis..."
	@cd data/ ; \
		./matlab do_analysis.m
	@mv ./data/Experiment_1.png ./docs/
	@mv ./data/Experiment_2.png ./docs/


document: ./docs/Experiment_1.png ./docs/Experiment_2.png 
	@echo "Knitting document"
	@R -e "xfun::pkg_attach2('rmarkdown')"
	@R -e "rmarkdown::render('docs/bae-LJC.Rmd')"
	@rm ./docs/*.log

