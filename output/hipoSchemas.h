#ifndef HIPO_SCHEMA_H
#define HIPO_SCHEMA_H 1


class HipoSchema {
public:
	HipoSchema();
	
	hipo::schema runConfigSchema;
	hipo::schema runRFSchema;
	hipo::schema trueInfoSchema;
	
	// generators
	hipo::schema geantParticle;
	hipo::schema mcEventHeader;
	hipo::schema userLund;
	hipo::schema lundParticle;

	// flux
	hipo::schema fluxADCSchema;

	// detectors
	hipo::schema alertAhdcADCchema;
	hipo::schema alertAtofADCchema;
	hipo::schema bandADCSchema;
	hipo::schema bandTDCSchema;
	hipo::schema bmtADCSchema;
	hipo::schema bstADCSchema;
	hipo::schema cndADCSchema;
	hipo::schema cndTDCSchema;
	hipo::schema ctofADCSchema;
	hipo::schema ctofTDCSchema;
	hipo::schema dcTDCSchema;
	hipo::schema dcDOCASchema;
	hipo::schema ecalADCSchema;
	hipo::schema ecalTDCSchema;
	hipo::schema fmtADCSchema;
	hipo::schema ftcalADCSchema;
	hipo::schema fthodoADCSchema;
	hipo::schema ftofADCSchema;
	hipo::schema ftofTDCSchema;
	hipo::schema ftrkTDCSchema;
	hipo::schema htccADCSchema;
	hipo::schema htccTDCSchema;
	hipo::schema ltccADCSchema;
	hipo::schema ltccTDCSchema;
	hipo::schema rfADCSchema;
	hipo::schema rfTDCSchema;
	hipo::schema richTDCSchema;
	hipo::schema rtpcADCSchema;
	hipo::schema rtpcPOSSchema;
	hipo::schema helADCSchema;
	hipo::schema helFLIPSchema;
	hipo::schema helONLINESchema;
	hipo::schema urwellADCSchema;
	hipo::schema rasterADCSchema;
	hipo::schema rawADCSchema;
	hipo::schema rawTDCSchema;
	hipo::schema rawSCALERSchema;
	hipo::schema rawVTPSchema;
	hipo::schema rawEPICSSchema;
	hipo::schema emptySchema;
	
	map<string, hipo::schema> schemasToLoad;
	
	// type: 0 = adc, 1 = tdc
	hipo::schema getSchema(string schemaName, int type) ;
	
	
};


#endif
