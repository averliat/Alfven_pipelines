#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup matplotlib for PDF output using PGF backend, add missing preambles from mplstyle.
"""

import matplotlib as mpl
#from matplotlib.backends.backend_pgf import FigureCanvasPgf
#mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
import matplotlib.pyplot as plt

plt.style.use('pdf')
mpl.rcParams['text.latex.preamble'] = [
        #r"\usepackage{nicefrac}",
        #r"\usepackage[lining,proportional]{ebgaramond}",
        #r"\usepackage{euler}",
        #r"\usepackage{siunitx}"
    ]
mpl.rcParams['pgf.preamble'] = [
        #r"\usepackage{nicefrac}",
        r"\usepackage{fontspec}",
        #r"\defaultfontfeatures[EBGaramond]{Ligatures=TeX,Numbers={Lining,Monospaced,Proportional}}",
        #r"\setmainfont[",
        #r"Extension   = .otf,",
        #r"UprightFont = *-Regular,",
        #r"ItalicFont  = *-Italic,",
        #r"BoldFont    = *-Bold.otf,",
        #r"BoldItalicFont = *-BoldItalic.otf",
        #r"]{EBGaramond}",
        #r"\usepackage[math-style=upright]{unicode-math}",
        #r"\setmathfont[Scale=MatchLowercase]{euler.otf}",
        #r"\setmathfont[Scale=MatchLowercase,range=\star]{latinmodern-math.otf}",
        #r"\usepackage{siunitx}"
    ]
