# Data

Several levels of data can be distinguished (see Figure below)
- Level 1: raw data in MAT format
- Level 2: raw data reorganized such that there is 1 file per fiber and 1 file containing force-pCa relations of all fibers (`force_pCa.mat`, provided in this folder). 
- Level 3: normalized data in which forces are normalized with respect to the force-pCa data
- Level 4: short-range stiffness (SRS) data

At this time, Level 3 data is provided to the manuscript reviewers and Level 4 data is included in this repository. 
The script `calc_data_SRS.m` convert Level 3 data into Level 4 data. Both the Level 3 data (containing time-series) and the Level 4 data (containing SRS estimates) can be used to reproduce the manuscript figures (see [biophysical-muscle-model/Reproduce](https://github.com/timvanderzee/biophysical-muscle-model/tree/main/Reproduce)).

[![](https://mermaid.ink/img/pako:eNp1U1uPojAU_ivNeTRIWgEZedhERd_mZWafdpmQI1QlgZbUMjf1v29bZ2QlI2lIe77L-U6THqGQJYcEtrV8K_aoNPmdZoIQMv87Gj1h05ISNZIxUfg2GhEWpLyYUBZj3i4xn7aa-g3qH-rM1X3f77HNHc3mRvNCXP_FoH-rZMEPB17eplBcqh2K6pOXQ9MB9OJsl8Z2XujqFXUlxd3h8ErJLWVo_QP8FTs1_mupCj42Mzl7Y2y-ra3ZOfsoq8GEQqoGaxv3JksvWA8Ez0_Phmn-fYhMOCaZk_H416m_gS_GiSwu-NLhBdZF_l-yE0m_DVJHuEa66lcXeHEfvhBWfQOL5SalxdfgwU5VJSRaddyDhhsHe4Sj1WWg97zhGSRmW_ItdrXOIBNnI2tR_JGy-VYq2e32kGyxPphT15ouPK1wp7C5VhUXJVdL2QkNCWPMmUByhHdIIurHIY2mcTijYciCqQcfkEzoxKfxNKbsIQomAYuCswefri31o2hmoDCk4UNklwe8rLRUj5cn5F7S-R9lmAhF?type=png)](https://mermaid.live/edit#pako:eNp1U1uPojAU_ivNeTRIWgEZedhERd_mZWafdpmQI1QlgZbUMjf1v29bZ2QlI2lIe77L-U6THqGQJYcEtrV8K_aoNPmdZoIQMv87Gj1h05ISNZIxUfg2GhEWpLyYUBZj3i4xn7aa-g3qH-rM1X3f77HNHc3mRvNCXP_FoH-rZMEPB17eplBcqh2K6pOXQ9MB9OJsl8Z2XujqFXUlxd3h8ErJLWVo_QP8FTs1_mupCj42Mzl7Y2y-ra3ZOfsoq8GEQqoGaxv3JksvWA8Ez0_Phmn-fYhMOCaZk_H416m_gS_GiSwu-NLhBdZF_l-yE0m_DVJHuEa66lcXeHEfvhBWfQOL5SalxdfgwU5VJSRaddyDhhsHe4Sj1WWg97zhGSRmW_ItdrXOIBNnI2tR_JGy-VYq2e32kGyxPphT15ouPK1wp7C5VhUXJVdL2QkNCWPMmUByhHdIIurHIY2mcTijYciCqQcfkEzoxKfxNKbsIQomAYuCswefri31o2hmoDCk4UNklwe8rLRUj5cn5F7S-R9lmAhF)

Note: Level 3 data is sufficient to reproduce the data portion of the manuscript. Optionally, `process_data.m` can be used to compute Level 2 and 3 data from Level 1 data, if the latter is provided. We intend to publish all data open-access upon manuscript publication. 

