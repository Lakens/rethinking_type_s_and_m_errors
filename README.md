# Rethinking Type S and M Errors

This repository contains files to computationally reproduce the manuscript "Rethinking Type S and M Errors" by Daniel Lakens, Cristian Mesquida, Gabriela Xavier-Quintais, Sajedeh Rasti, Enrico Toffalini, and Gianmarco Altoé. For questions, email D.Lakens\@tue.nl


# Abstract

Gelman and Carlin (2014) introduced Type S (sign) and Type M (magnitude) errors to highlight the possibility that statistically significant results in published articles are misleading. While these concepts have been proposed for both when designing a study (prospective) and when evaluating results in a published paper (retroactive), we argue that these statistics do not facilitate the proper design of studies, nor the meaningful interpretation of results. Type S errors are a response to the criticism of testing against a point null of exactly zero in contexts where true zero effects are implausible. Testing against a minimum effect size, while controlling the Type 1 error rate, provides a more coherent and practically useful alternative. Type M errors are primarily a tool to warn against effect size inflation after selectively reporting significant results, but we argue that statistical indices such as the critical effect size or bias adjusted effect size are preferable approaches. We do believe that Type S and M errors can be valuable pedagogically as part of statistics education where the principles of error control are explained, and consequences of bad research practices are explained. Type S and M errors could be used in the discussion section of studies that fail to follow good research practices. Overall, we argue their use-cases are more limited than is currently recognized, and alternative solutions deserve greater attention.

# Repository content: 

- rethinking_type_s_and_m_errors.qmd (Quarto document containing the reproducible manuscript)
- rethinking_type_s_and_m_errors.pdf (the manuscript)
- session_info.txt (text file describing all packages used when compiling the manuscript)
- n=40d=0.5sim=1000.rds (the results of the simulation, saved as R data file to save time when compiling the manuscript)
- references.bib (a bib file containing all references)
- README.md (this readme file)
- An _extensions folder related to the apaquarto extension
- A folder rethinking_type_s_and_m_errors_files generated while rendering the manuscript
- A license file


# Compiling

Make sure the working directory is set to the rethinking_type_s_and_m_errors.qmd folders. 
Typst is used to create the PDF file. I have changed the font in "rethinking_type_s_and_m_errors\_extensions\wjschne\apaquarto\typst\typst-template.typ" on line 30, if you do not have Segoe UI installed, you might need to change it. 