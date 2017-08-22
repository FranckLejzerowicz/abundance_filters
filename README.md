# abundance_filters

Mutlidimension abundance filter: per sample, across dataset, on selected or fractions of samples or replicates.

This script offers various options that can be combined to filter an OTU table using an abundance threshold. The filtering applies for each OTU independantly but it may be set to filter the OTU according to various features of its distribution in the samples. The most simple filter will remove the OTU presence in each sample independently ('-meth simple'), or across the entire dataset if the option '--sum' is added. When the option '--sum' is activated, it is the sum of the abundances of an OTU in all selected samples that is compared to the threshold, and not the abundance of the OTU in each sample (or group of samples). These "selected samples" are a subset of the input table samples, which can be defined as a fraction / a percent of the samples, or a selection using indexing and/or regular expressions, or specified replicates. The option '--only' will make the filter apply only on the "selected samples" (e.g. if the OTU is abundant enough in the samples selected using '-meth choice'), and will return zero count for the other samples. It is recommended to use '--only' in conjonction with the option '--selection', that will write the output table with only the "selected samples" in columns (i.e. avoids writing out zero-filled columns).

It is mandatory that the diversity-to-sample input table has samples in columns and diversity entity entries in rows. Each row must contain distribution data (e.g. numbers of sequences) for an entry in the samples. An entry can be an OTU, an ISU (Individual Sequence Unit or unique sequences), a taxon, or any other similar entity that makes up the table.

The table passed as input to the '-i' option is automatically checked to detect (i) if the first line is a header containing the samples names, (ii) the number of left-hand columns that do not correspond to numeric data and hence will be interpreted as entry metadata (e.g. OTU name, taxonomy, OTU pH, ...). Alternatively, whether (i) is true can be set using option '--head' and the number (ii) of fields can be set using option '-meta'. In the case where the table has no left-hand metadata (i.e. only numbers), a name can be created for each OTU if a prefix is entered using option '-name_otus'.


## Usage
```
usage: abundance_filter.py [-h] -i I [-o [O]] [-t [T]] [-x [X [X ...]]] [-meth [{simple,minimum,presence,choice,replicates}]] [-mode [MODE [MODE ...]]] [--sum] [--only] [--selection] [-name-otus NAME_OTUS] [-f [F]] [-sep [SEP]] [--h]
```

### Optional arguments
```
  -h, --help            show help message and exit
  -i I                  Input OTU/ISU table name (required)
  -o [O]                Output OTU/ISU table name (default = add arguments extensions to input name, and '.stats' for stats file)
  -t [T]                Abundance threshold, inclusive (default = 10)
  -x [X [X ...]]        Samples to apply filter on and to output, using:
                            * indices numbers (from 1 to number of samples)
                            * regex for sample(s) names)
  -meth [{simple,minimum,presence,choice,replicates}]
                        Filtering method: [default = 'simple']
                            * 'simple': filter OTUs/ISUs not abundant enough
                            * 'minimum': filter OTU/ISU not abundant enough in a percent of samples. Use with '-mode' [default = 10]
                            * 'presence': filter OTUs/ISUs not present in a min percent of samples. Use with '-mode' [default = 10]
                            * 'choice': define groups of samples where OTUs/ISUs must be abundant enough. Use with -mode
                            * 'replicates': use groups of samples among which OTUs/ISUs must be abundant enough. Use with -mode)
                            [EACH METHOD CAN OPERATE PER SAMPLE OR ACROSS THE ENTIRE DATASET, USING '--sum']
  -mode [MODE [MODE ...]]
                        Filtering mode, for '-meth':
                            * 'minimum': percent of samples (integer between 1 and 100) [default = 10]
                        	(remove OTU/ISU that do not reach threshold in the percent of samples)
                        	e.g. '-meth minimum -mode 50' = OTU/ISU must be abundant enough in at least half of the samples
                            * 'presence': percent of samples (integer between 1 and 100) [default = 10]
                        	(remove OTU/ISU that is not present in the percent of samples)
                        	e.g. '-meth presence -mode 50' = OTU/ISU must be present in at least half of the samples
                           (** Both '-meth choice' and '-meth replicates' could be combined with a minimum percent of these samples [default = 100])
                            * 'choice':
                        	- one (or more) space-separated regex to select by sample name
                        	- at least 2 comma-separated integers to select by sample indices (starts at 1)
                        	- both regex and indices
                        	[See '--select' for output]
                        	e.g. '-meth choice -mode gut_.*_100ppm 30' = OTU/ISU must be abundant enough in 30 % of the regex-matching samples
                            * 'replicates': file with tab-/comma-separated replicates samples names in rows
                        	Optional integers for:
                        	- number of replicate occurrences per group of sample replicates [default = 2]
                        	- number of occurrences in terms of groups of samples replicates (or percent if float between 0 and 1) [default = ]
                        	e.g. '-meth replicates -mode <repFile.txt> 3 10' = OTU/ISU must be abundant enough in >=3 replicates of >= 10 groups of replicates
                           ['minimum', 'choice' and 'replicates': see option '--only' to apply OTU filtering in selected samples only]
  --sum                 Compare to the threshold the SUM of the abundances in all the samples to be filtered.
                        (default = not activated = filter per sample)
  --only                Remove ISUs/OTUs only in the samples where they do not reach the abundance threshold.
                        Automatic with '-meth sample' / Manual for methods 'minimum', 'choice' and 'replicates'.
                        (default = remove OTU in all samples, or in all replicates of a group of replicates for method 'replicates').
  --selection           Return filtered data table containing only the samples selected using '-meth choice'.
  -name-otus NAME_OTUS  Character string for pass-filter OTUs/ISUs prefix, with OTUs/ISUs indices as suffix
                        (only if input table contains no name, else default = no name if option not used).
  -meta [META]          Number of non-numeric column in the OTU/ISU table (default = Auto-detect - Warning if OTU/ISU names are numeric)
  --head                First line is a header line (default = automatic detection).
  -sep [SEP]            Fields separator for the input table (default = tabulation).
```


#### Requirements
Python2.7<br />
Numpy<br />
Biopython<br />
R<br />

