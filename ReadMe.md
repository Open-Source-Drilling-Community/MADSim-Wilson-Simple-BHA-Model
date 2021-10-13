MADSim-Wilson-Simple-BHA-Model: Simplified Nonlinear BHA Model

The supplied code is written in Matlab (versions 2013 - 2021). It is based on a fully nonlinear downhole dynamics model, and 
has been simplified to provide insight to the implementation of nonlinear finite elements as they relate to BHA mechanics. A
general overview of the model is provided in the document "Overview_MADSim, SDI OpenSource.pdf".

Provided files:

SDIopneSource_MADSim.m						Source code, to run the model use this file in MATLAB (no additional toolboxes should be required)
Overview_MADSim, SDI OpenSource.pdf				General overview of provided model, references contained within this document
LICENSE								License, modeled after MIT open-source license
BHA Data_Example RSS.xlsx					Example data, see "Overview_MADSim, SDI OpenSource.pdf"
BHA Data_ExampleAssembly.xlsx					Example data, see "Overview_MADSim, SDI OpenSource.pdf"
Example Results_Loading_ExampleAssemblyTurn.xlsx		Example data, see "Overview_MADSim, SDI OpenSource.pdf"
Example Results_MotorYield (Sliding)_Example Assembly.xlsx	Example data, see "Overview_MADSim, SDI OpenSource.pdf"
Wilson, 2017 (Dissertation).pdf					Reference to original model, from which this simplifeid model is based


Limitations of provided model:

+	Straight wellbores only
+	Constant wellbore OD
+	No buckling*
+	No dynamic calculations
+	Graphic user interfaces have been removed
+	Advanced plotting features have been removed



Due to copyright limitations, only some of the references listed for this Open-Source contribution can be provided here.

1)	Wilson, J.K. 2017. “Nonlinear Drillstring Modeling with Application to Induced Vibrations in Unconventional Horizontal Wells”. Ph.D. Dissertation. Texas A&M University.
	
	Updated document provided. Link to original document: https://oaktrust.library.tamu.edu/handle/1969.1/161300

2)	Wilson, J.K. 2019. “Field Validation of a New Bottomhole-Assembly Model for Unconventional Shale Plays”. SPE 191780 PA. SPE Drilling & Completion.

	Link to document: https://onepetro.org/DC/article/34/02/189/207025/Field-Validation-of-a-New-Bottomhole-Assembly?searchresult=1

3)	Heisig, G. 1995. “Postbuckling Analysis of Drillstrings Using the Finite-Element Method”. ASME Drilling Technology, Vol.65.

	Link to document record: https://www.osti.gov/biblio/170144
	Document unavailable online. Access can usually be gained through a University library

4) 	Wilson, J.K. 2018. "Field Validation of a New BHA Model and Practical Case Studies in Unconventional Shale Plays, with a Framework for Automated Analysis for Operations Support". SPE-191780-18ERM-MS

	Link to document: https://onepetro.org/SPEERM/proceedings/18ERM/3-18ERM/D033S001R004/449253?searchresult=1
