# Data

Several levels of data can be distinguished (see Figure below)
- Level 1: raw data in MAT format
- Level 2: raw data reorganized such that there is 1 file per fiber and 1 file containing force-pCa relations of all fibers (`force_pCa.mat`, provided in this folder). 
- Level 3: normalized data in which forces are normalized with respect to the force-pCa data

At this time, Level 3 data is provided to the manuscript reviewers. 
This data can be used to compute short-range stiffness (summary) data using `calc_data_SRS.m`.
The resulting file is called `SRS_data.mat` and is provided in this folder. 
Both the Level 3 data (containing time-series) and the SRS summary data can be used to reproduce the manuscript figures (see [biophysical-muscle-model/Reproduce](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce)).

[![](https://mermaid.ink/img/pako:eNp1U0tvgkAQ_iubORokuwiIHJqo2FsvTU8thoywKgmwZF1q6-O_d1lflSgnZr75XgnsIRUZhxCWhdima5SKfERxRQgZf_V671jWJEOFpE8kbns9wgYRTx3KhpjUU0z8WlG7RPVgz8zetu0btnjCWdxx5sT4Tzr-tRQp32x4dp9CciFXWOU7nnVFO9DcyE617DhV-TeqXFRPy-H1JGlPutIP4HPsSOu_Cpnyvu5k5LWwfpbtru15izLrNKyELLFo495lORMMhYxJv_9yuDU72x_I5IRPDZ5ikSb_HA8kughE5uBqdeXPTvDkGQwWrGSeQahkwy0oub5oR9i3xBjUmpc8hlC_ZnyJTaFiiKujptVYfQpRXphSNKs1hEssNnpqau3AoxxXEsvrVvIq43IqmkpB6LrUiEC4hx8IGaP2iPp-MPQDx_EC17fgF0InYDbV02jAqOeOmHe0YGdsqT0IfM8fjJjjDv0h9RwLeJYrId9On775A45_VNnyWw?type=png)](https://mermaid.live/edit#pako:eNp1U0tvgkAQ_iubORokuwiIHJqo2FsvTU8thoywKgmwZF1q6-O_d1lflSgnZr75XgnsIRUZhxCWhdima5SKfERxRQgZf_V671jWJEOFpE8kbns9wgYRTx3KhpjUU0z8WlG7RPVgz8zetu0btnjCWdxx5sT4Tzr-tRQp32x4dp9CciFXWOU7nnVFO9DcyE617DhV-TeqXFRPy-H1JGlPutIP4HPsSOu_Cpnyvu5k5LWwfpbtru15izLrNKyELLFo495lORMMhYxJv_9yuDU72x_I5IRPDZ5ikSb_HA8kughE5uBqdeXPTvDkGQwWrGSeQahkwy0oub5oR9i3xBjUmpc8hlC_ZnyJTaFiiKujptVYfQpRXphSNKs1hEssNnpqau3AoxxXEsvrVvIq43IqmkpB6LrUiEC4hx8IGaP2iPp-MPQDx_EC17fgF0InYDbV02jAqOeOmHe0YGdsqT0IfM8fjJjjDv0h9RwLeJYrId9On775A45_VNnyWw)

Note: Level 3 data is sufficient to reproduce the data portion of the manuscript. Optionally, `process_data.m` can be used to compute Level 2 and 3 data from Level 1 data, if the latter is provided. We intend to publish all data levels open-access upon manuscript publication. 

