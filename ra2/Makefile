# Make eps and pdf from combined eps/LaTeX xfig graphics


all: js_DeltaPhi_vs_MET.eps js_Smearing.eps

js_DeltaPhi_vs_MET.eps: js_DeltaPhi_vs_MET.pstex js_DeltaPhi_vs_MET.pstex_t
	latex js_DeltaPhi_vs_MET.tex
	dvips js_DeltaPhi_vs_MET.dvi
	ps2epsi js_DeltaPhi_vs_MET.ps
	mv js_DeltaPhi_vs_MET.epsi js_DeltaPhi_vs_MET.eps
	epstopdf js_DeltaPhi_vs_MET.eps
	rm js_DeltaPhi_vs_MET.aux js_DeltaPhi_vs_MET.log js_DeltaPhi_vs_MET.dvi js_DeltaPhi_vs_MET.ps

js_Smearing.eps: js_Smearing.pstex js_Smearing.pstex_t
	latex js_Smearing.tex
	dvips js_Smearing.dvi
	ps2epsi js_Smearing.ps
	mv js_Smearing.epsi js_Smearing.eps
	epstopdf js_Smearing.eps
	rm js_Smearing.aux js_Smearing.log js_Smearing.dvi js_Smearing.ps
