%% This file is 'aes2e.cls' 
%%
%% Created by Aptara Inc.
%% Dated 26 October 2010
%%
%% (c) 2010 Audio Engineering Society (AES)
%%
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version 1.0 (History)
%% Released on 26 October 2010
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{aes2e}[2010/10/16 v1.0 Standard LaTeX document class for AES]
%
\newif\ifCUPmtlplainloaded
%
\@twosidetrue\@mparswitchtrue\@twocolumnfalse
%
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\ExecuteOptions{final}
\ProcessOptions
%
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\bibname{REFERENCES}
\newcommand\indexname{Index}
\newcommand\figurename{Fig.}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\chaptername{Chapter}
\newcommand\appendixname{Appendix}
%
\newdimen\trimwidth
\newdimen\trimheight
\newdimen\typewidth
\newdimen\typeheight
%
\newdimen\tempdimen
\newdimen\normaltextheight
\newskip\normalbaselineskip
\newdimen\blindfoliodrop
%
\newdimen\enumdim
\newdimen\figwidth
\newdimen\tabledim
\newskip\tableleftskip
\newskip\tablerightskip
\newskip\tabnoteleftskip
\newskip\tabnoterightskip
%
\newdimen\mathindent
%
\newbox\tempbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Boxes & Dimensions %%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\RequirePackage{mathptmx}
\RequirePackage{helvet}
\RequirePackage{courier}
% headings
\def\rhfont{\fontfamily{\sfdefault}\fontsize{7}{8}\selectfont}
\def\RHfont{\fontfamily{\sfdefault}\fontsize{7}{8}\selectfont}
\def\foliofont{\fontfamily{\sfdefault}\fontsize{7}{8}\selectfont}
\def\dropfoliofont{\fontfamily{\sfdefault}\fontsize{7}{8}\selectfont}
\def\sptitlefont{\fontsize{14}{15}\rightskip3pc plus1fill\selectfont\bfseries}
\def\titlefont{\fontfamily{\sfdefault}\fontsize{20}{24}\leftskip0pt plus1fill\rightskip\leftskip\selectfont\bfseries}
\def\authorfont{\reset@font\fontfamily{\sfdefault}\fontsize{10}{10}\selectfont\bfseries\centering}
\def\authorrolefont{\reset@font\fontfamily{\sfdefault}\fontsize{10}{10}\selectfont\itshape}
\def\authoremailfont{\reset@font\fontsize{8}{8}\selectfont}%
\def\affilfont{\reset@font\fontsize{10}{12}\rightskip0pc plus1fill\leftskip\rightskip\selectfont\itshape}
\def\abstractfont{\fontsize{9}{10}\leftskip6.25pc\rightskip\leftskip\selectfont}
\def\keywordsfont{\rightskip0pt plus1fill}
%% section
\def\sectionfont{\fontfamily{\sfdefault}\fontsize{10}{12}\rightskip0pc plus1fill\selectfont\bfseries}
\def\subsectionfont{\fontfamily{\sfdefault}\fontsize{10}{12}\rightskip0pc plus1fill\selectfont\bfseries}
\def\subsubsectionfont{\fontfamily{\sfdefault}\fontsize{10}{12}\rightskip0pc plus1fill\selectfont\bfseries\itshape}
\def\paragraphfont{\itshape}%
\def\subparagraphfont{}%
\def\subsubparagraphfont{}%
%% figure
\def\figcaptionfont{\fontsize{9}{10}\selectfont}
\def\figcaptionnumfont{\fontsize{9}{10}\selectfont}
\def\tablefont{\fontsize{9}{10}\selectfont}%
\def\tablecaptionfont{\fontsize{9}{10}\selectfont\leftskip\tableleftskip\rightskip\tablerightskip}% 
\def\tablecaptionnumfont{\fontsize{9}{10}\selectfont}% 
\def\tabnotefont{\fontsize{8}{9}\selectfont\leftskip\tabnoteleftskip\rightskip\tabnoterightskip}%
%% box
\def\boxfont{\fontsize{8}{9}\selectfont}%
\def\boxcaptionfont{\fontsize{8}{9}\selectfont}% 
\def\boxcaptionnumfont{\fontsize{8}{9}\selectfont\bfseries}% 
%% BM
\def\indexfont{\fontsize{9}{11}\selectfont\raggedright}
\def\bibliofont{\fontsize{10}{12}\selectfont}
%
\def\listfont{}
\def\listdevicefont{}
\def\extractfont{\fontsize{9}{11}\selectfont\leftskip1pc\rightskip\leftskip}
%
%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
\normalbaselineskip12pt
%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xpt{12}
   \abovedisplayskip 5.5\p@ \@plus2\p@ \@minus.5\p@
   \abovedisplayshortskip \abovedisplayskip
   \belowdisplayshortskip \abovedisplayskip
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\def\biggg#1{{\hbox{$\left#1\vbox to20.5\p@{}\right.\n@space$}}}
\def\Biggg#1{{\hbox{$\left#1\vbox to23.5\p@{}\right.\n@space$}}}
\normalsize
%
\newcommand\small{%
   \@setfontsize\small\@ixpt{11}%
   \abovedisplayskip 8.5\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setSmallDelims}
%
\def\setSmallDelims{%
\def\big##1{{\hbox{$\left##1\vbox to7.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to10.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to13.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to16.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to19.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to22.5\p@{}\right.\n@space$}}}%
}
%
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@ixpt{11}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
   \setFootnotesizeDelims}
%
\def\setFootnotesizeDelims{%
\def\big##1{{\hbox{$\left##1\vbox to6.5\p@{}\right.\n@space$}}}%
\def\Big##1{{\hbox{$\left##1\vbox to9.5\p@{}\right.\n@space$}}}%
\def\bigg##1{{\hbox{$\left##1\vbox to12.5\p@{}\right.\n@space$}}}%
\def\Bigg##1{{\hbox{$\left##1\vbox to15.5\p@{}\right.\n@space$}}}%
\def\biggg##1{{\hbox{$\left##1\vbox to18.5\p@{}\right.\n@space$}}}%
\def\Biggg##1{{\hbox{$\left##1\vbox to21.5\p@{}\right.\n@space$}}}%
}
%   
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
%
\newcommand{\cal}{\mathcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
\newcommand{\mit}{\protect\pmit}
\newcommand{\pmit}{\@fontswitch{\relax}{\mathnormal}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\trimheight{10.875in}
\setlength\trimwidth{8.0625in}

\setlength\paperheight {10.875in}%
\setlength\paperwidth  {8.0625in}%
\def\special@paper{8.0625in,10.875in}
\special{papersize=8.0625in,10.875in}
%
\setlength\typeheight{61\baselineskip}
\setlength\typewidth{40.5pc}
%
\setlength\textheight{59\baselineskip}
\addtolength\textheight{-4.15pt}
\setlength\normaltextheight{\textheight}
\setlength\textwidth{40.5pc}
%
\setlength\topmargin{22.5pt}
\setlength\oddsidemargin   {54pt}
\setlength\evensidemargin\trimwidth
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
%
\setlength\headheight{6\p@}% 
\setlength\headsep   {14.6pt}% \typeheight - \textheight - \headheight
\setlength\topskip   {7.6\p@}
%
\setlength\footskip{21.5pt}%.35in}
\setlength\maxdepth{4\p@}
%
\blindfoliodrop\trimheight
\advance\blindfoliodrop-\typeheight
\advance\blindfoliodrop-\topmargin
\advance\blindfoliodrop-\footskip
\advance\blindfoliodrop18pt
%
\setlength\parindent{1em}
%
\setlength\marginparwidth {1pc}%.75in}
\setlength\marginparsep{3\p@}
\setlength\marginparpush{5\p@}
%
\setlength\footnotesep{6.65\p@}
\setlength{\skip\footins}{21\p@ \@plus 8\p@ \@minus 8\p@}% change + based on FN font height
\skip\@mpfootins = \skip\footins
%
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@}
\setlength\@fpsep{8\p@ \@plus 1fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@}
\setlength\@dblfpsep{8\p@ \@plus 1fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
%
\setlength\partopsep{0pt}
\setlength\lineskip{1\p@}% check if it can be flexible
\setlength\normallineskip{1\p@}%
\renewcommand\baselinestretch{}
\setlength\parskip{\z@}%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Parameters %%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Article Title  %%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\sptitle#1{\gdef\@sptitle{#1}}
\sptitle{}
\def\title#1{\gdef\@title{#1}}
\def\authorgroup#1{\long\gdef\@author{\par\noindent{\authorfont#1\par}}}%
\def\author#1{{\authorfont#1}}
\def\role#1{{\authorrolefont#1}}
\def\email#1{{\vskip1.8pt\authoremailfont#1\par}}
\def\affil#1{\par\vskip12pt\reset@font\affilfont#1\vphantom{pqj}\par\vskip12pt}
%
\newif\ifauthorgd
\renewcommand\thefootnote{\@fnsymbol\c@footnote}%
\newcommand\maketitle{\clearpage\par%
  \begingroup%
    \ifauthorgd\@maketitle\else
    \twocolumn[\@maketitle]\fi%
    \thispagestyle{plain}%
  \@thanks%
  \endgroup%
  \let\@maketitle\relax
   \gdef\@author{}%
  \gdef\@title{}%
  \gdef\@sptitle{}%
  \gdef\@thanks{}%
  \renewcommand\thefootnote{\arabic{footnote}}%
  \setcounter{footnote}{0}%
  \setcounter{section}{-1}%
}
%
\def\@maketitle{%
  \clearpage%
  \null%
  \bgroup%
  \vspace*{8pt}%
    {\ifx\@title\@empty\else\titlefont\@title\vphantom{pg}\fi\par\vspace*{35pt}}%
    \@author\vskip26pt
    \setlength{\parindent}{1em}
    {\ifx\@abstract\empty\else{\reset@font\abstractfont\@abstract\par}\fi}%
     \vskip50.5pt%
  \egroup}
%
%=========================================================
%
\def\qdata#1#2{\par\noindent\@hangfrom{\hbox to 5.5pc{#1\hfill}}{#2\par\vskip\baselineskip}}
%
\newenvironment{conf}{\par\addvspace{11pt}%
                \parindent0pt\ignorespaces%
               }{\par\addvspace{0pt}}
%
\def\@abstract{}%
\def\abstract#1{\gdef\@abstract{\par#1\par}}%
%
\def\keyword#1{#1}
\newenvironment{keywords}{\par\addvspace{11pt}
                \noindent\keywordsfont\ignorespaces%
               }{\par\addvspace{0pt}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sectioning commands  %%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{4}
\newcounter {part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart          {\Roman{part}}
\renewcommand\thechapter       {\arabic{chapter}}
\renewcommand\thesection       {\arabic{section}}
\renewcommand\thesubsection    {\thesection.\arabic{subsection}}
\renewcommand\thesubsubsection {\thesubsection .\arabic{subsubsection}}
\renewcommand\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph  {\theparagraph.\arabic{subparagraph}}
%
\def\@seccntformat#1{\csname the#1\endcsname\hskip5pt}
%
\newcommand\section{\@startsection {section}{1}{\z@}{16pt \@plus 2pt \@minus .2pt}{8pt}{\sectionfont}}%
\newcommand\subsection{\@startsection{subsection}{2}{\z@}{12pt \@plus 2pt \@minus .2pt}{2pt}{\subsectionfont}}%
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{12pt \@plus 2pt \@minus .2pt}{2pt}{\subsubsectionfont}}%
\newcommand\paragraph{\@startsection{paragraph}{4}{1em}{-0.0001pt}{-1em}{\paragraphfont}}%
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}{-19pt \@plus -2pt \@minus -.2pt}{-1em}{\subparagraphfont}}%
\newcommand\subsubparagraph{\@startsection{subsubparagraph}{6}{\parindent}{-19pt \@plus -2pt \@minus -.2pt}{-1em}{\subsubparagraphfont}}%
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth%
    \let\@svsec\@empty%
  \else%
    \refstepcounter{#1}%
    \let\@@protect\protect
    \def\protect{\noexpand\protect\noexpand}%
    \edef\@svsec{\@seccntformat{#1}}%
    \let\protect\@@protect%
  \fi
  \@tempskipa #5\relax%
  \ifdim \@tempskipa>\z@%
    \begingroup #6\relax%
      \ifx\disp@level\DISPLVL%
        \noindent{\hskip #3\relax\@svsec}%
              {\interlinepenalty \@M {#8\hfill}\par}%
      \else%
         \noindent{\hskip #3\relax\@svsec}%
              {\interlinepenalty \@M #8\par}%
     \fi%
    \endgroup%
   \ifx\sec@@type\empty\csname #1mark\endcsname{#7}\addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else<numberline>\@seclabeltoc</numberline>\fi%
      #7}\fi%
  \else%
    \def\@svsechd{#6\hskip #3\relax%
      \@svsec #8\ifnum#2=3\else\ifnum#2=4\fi\fi\csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else\protect\numberline{\csname the#1\endcsname}\fi#7}}%
  \fi%
  \@xsect{#5}%
\gdef\@sectionlabel{}}%
%
\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec \leavevmode \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim \@tempskipa <\z@
    \@tempskipa -\@tempskipa \@afterindentfalse
  \fi
  \if@nobreak
  \ifnum#2=3
  \fi
    \everypar{}
  \else
    \addpenalty\@secpenalty\addvspace\@tempskipa
  \fi
  \@ifstar
    {\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%  \@ssect{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}{ARG}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Sectioning commands  %%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\topsepi \topsepi6\p@ \@plus2\p@%
\newskip\topsepii \topsepii2pt%
\newskip\topsepiii \topsepiii2pt%
\newskip\itemsepi \itemsepi0pt
\newskip\itemsepii \itemsepii0pt
\newskip\itemsepiii \itemsepiii0pt
\newdimen\LabelSep \LabelSep4.7pt
%
\def\@listI{\leftmargin\leftmargini
            \labelwidth\leftmargini
            \advance\labelwidth-\labelsep
            \parsep 0\p@%
            \topsep \topsepi
            \itemsep\itemsepi}%
\let\@listi\@listI
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep\topsepii
              \parsep 0pt
              \itemsep\itemsepii}
\def\@listiii {\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep\topsepiii
              \parsep 0pt
              \itemsep\itemsepiii}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\setlength\leftmargini  {2.5em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\setlength\leftmargin    {\leftmargini}
%
\setlength  \labelsep  {\LabelSep}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
%
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\alph{enumii}}
\renewcommand\theenumiii{\roman{enumiii}}
\renewcommand\theenumiv{\Alph{enumiv}}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{\theenumii.}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\font\lcir = lcircle10 at 10pt%
\newcommand\bulls{\raise1.5pt\hbox{\lcir\char'162}}
\def\textbullet{\leavevmode\raise1.2pt\hbox{\bulls}}
\def\textendash{textbf{--}}
\def\textasteriskcentered{\leavevmode\raise-1.5pt\hbox{*}}
\def\textperiodcentered{\leavevmode\raise1.5pt\hbox{\bulls}}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
%
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
%
\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\itemsep\z@
                        \itemindent-1.5em%
                        \listparindent\itemindent
                        \rightmargin\leftmargin
                        \advance\leftmargin 1.5em}%
                \item[]}
               {\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Quotes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{quote}{\list{}{\rightmargin\leftmargin}%
                \item[]}
               {\endlist}
%
\newenvironment{extract}{\par\addvspace{13.5pt plus2pt}
      \noindent\extractfont\ignorespaces 
}{%
      \par\addvspace{12pt plus2pt}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Quotes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Enumerate list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newskip\listtextleftmargin\listtextleftmargin 0pt%
\newskip\listtextleftmarginii\listtextleftmarginii0pt%
\newskip\listtextleftmarginiii\listtextleftmarginiii0pt%
\newskip\listtextrightmargin\listtextrightmargin0pt%
\newskip\listlabelleftskip \listlabelleftskip0pt%
\newskip\listlabelleftskipii \listlabelleftskipii0pt%
\newskip\listlabelleftskipiii \listlabelleftskipiii0pt%
\newskip\abovelistskipi\abovelistskipi12pt%
\newskip\belowlistskipi\belowlistskipi12pt%
\newskip\abovelistskipii\abovelistskipii0pt%
\newskip\belowlistskipii\belowlistskipii0pt%
\newskip\abovelistskipiii\abovelistskipiii0pt%
\newskip\belowlistskipiii\belowlistskipiii0pt%
\newskip\labelsepi \labelsepi5pt
\newskip\labelsepii \labelsepii5.2pt
\newskip\labelsepiii \labelsepiii\z@
\newdimen\enumdimwd
\newif\iflabelrightalign\labelrightaligntrue
%
\def\listdevicefonti{}
\def\listdevicefontii{}
\def\listdevicefontiii{}
\def\listfont{}
\def\textlistlabel{}
%
\def\enummax#1{%
  \labelsep\csname labelsep\romannumeral\the\@enumdepth\endcsname
  \ifdim\listtextleftmargin>\z@\labelsepi0pt\fi
  \ifdim\listtextleftmarginii>\z@\labelsepii0pt\fi
  \ifdim\listtextleftmarginiii>\z@\labelsepiii0pt\fi
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1\hskip\labelsep}%
  \enumdim\wd\tempbox
  \setbox\tempbox\hbox{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname#1}%
  \enumdimwd\wd\tempbox
  \expandafter\global\csname leftmargin\romannumeral\the\@enumdepth\endcsname\enumdim
  \ifdim\listtextleftmargin>\z@
     \leftmargini\listtextleftmargin
     \ifdim\listlabelleftskip>\z@
        \advance\leftmargini-\listlabelleftskip
     \fi  
  \fi
  \ifdim\listtextleftmarginii>\z@
     \leftmarginii\listtextleftmarginii
     \ifdim\listlabelleftskipii>\z@
        \advance\leftmarginii-\listlabelleftskipii
     \fi  
  \fi
  \ifdim\listtextleftmarginiii>\z@
     \leftmarginiii\listtextleftmarginiii
     \ifdim\listlabelleftskipiii>\z@
        \advance\leftmarginiii-\listlabelleftskipiii
     \fi  
  \fi
  \ifdim\listlabelleftskip>\z@
     \advance\leftmargini\listlabelleftskip
  \fi
  \ifdim\listlabelleftskipii>\z@
     \advance\leftmarginii\listlabelleftskipii
  \fi
    \ifdim\listlabelleftskipiii>\z@
     \advance\leftmarginiii\listlabelleftskipiii
  \fi}
%
\enummax{1.}
%
\def\enumerate{\@ifnextchar[{\@enumerate}{\@enumerate[\csname label\@enumctr\endcsname]}}%% 
%
\def\@enumerate[#1]{\par
      \ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth\@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \setcounter{\@enumctr}{1}\enummax{#1}%
      \list
                                                    {\csname label\@enumctr\endcsname}{\usecounter{\@enumctr}%
      \topsep\csname abovelistskip\romannumeral\the\@enumdepth\endcsname
      \itemsep\csname itemsep\romannumeral\the\@enumdepth\endcsname
      \listfont %
      \ifnum \@enumdepth=1 
         \rightmargin\listtextrightmargin
         \advance\rightmargin\rightskip
         \advance\leftmargin\leftskip
         \tempdimen\leftmargini
         \advance\tempdimen-\labelsep
	 \ifparalist%
\def\makelabel##1{{\hspace*{\listlabelleftskip}{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname{\textlistlabel##1}}}}%%
\else%
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskip{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi\textlistlabel##1}}}}%
\fi%
      \else
      \ifnum \@enumdepth=2
         \tempdimen\leftmarginii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
      \else
         \ifnum \@enumdepth=3
         \tempdimen\leftmarginiii
         \advance\tempdimen-\labelsep
         \def\makelabel##1{\hbox to \tempdimen{\hskip\listlabelleftskipiii{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname\hbox to \enumdimwd{\iflabelrightalign\hss\fi##1}}}}%
     \else
         \def\makelabel##1{\hss\llap{\csname listdevicefont\romannumeral\the\@enumdepth\endcsname##1}}%
      \fi\fi
      \fi}\fi}
%
\def\endenumerate{\@topsepadd\csname belowlistskip\romannumeral\the\@enumdepth\endcsname\endlist}%
%
\newif\ifparalist\global\paralistfalse%
\newenvironment{paralist}{%
\global\paralisttrue%
 \labelsepi5.5pt%
\listtextleftmargin.001pt%
\listlabelleftskip10pt%  
  \def\theenumi{\arabic{enumi}}
  \def\theenumii{\arabic{enumii}}%
  \def\theenumiii{\arabic{enumiii}}%
  \def\theenumiv{\arabic{enumiv}}%
  \renewcommand\labelenumi{\theenumi)}
  \begin{enumerate}%
}{%
  \end{enumerate}\global\paralistfalse}
%%
\newenvironment{arabiclist}{%
\labelsepi5.5pt%
\listlabelleftskip10pt%  
  \def\theenumi{\arabic{enumi}}
  \def\theenumii{\arabic{enumii}}%
  \def\theenumiii{\arabic{enumiii}}%
  \def\theenumiv{\arabic{enumiv}}%
  \renewcommand\labelenumi{\theenumi)}
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{wherelist}{%
  \def\listdevicefonti{\upshape}
  \listtextleftmargin 0pt
  \abovelistskipi13\p@ plus2pt
  \belowlistskipi13\p@ plus2pt
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{romanlist}{%
  \def\theenumi{\roman{enumi}}\def\theenumii{\roman{enumii}}%
  \def\theenumiii{\roman{enumiii}}\def\theenumiv{\roman{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{alphalist}{%
\listtextleftmargin 18pt
\listtextleftmarginii 18pt
\itemsepi0pt
\labelrightalignfalse
  \def\theenumi{\alph{enumi}}\def\theenumii{\alph{enumii}}%
  \def\theenumiii{\alph{enumiii}}\def\theenumiv{\alph{enumiv}}%
  \renewcommand\labelenumi{\theenumi)}
  \renewcommand\labelenumii{\theenumii)}
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{Romanlist}{%
  \def\theenumi{\Roman{enumi}}\def\theenumii{\Roman{enumii}}%
  \def\theenumiii{\Roman{enumiii}}\def\theenumiv{\Roman{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{Alphalist}{%
  \def\theenumi{\Alph{enumi}}\def\theenumii{\Alph{enumii}}%
  \def\theenumiii{\Alph{enumiii}}\def\theenumiv{\Alph{enumiv}}%
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  unnumlist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{unnumlist}{\par\addvspace{12pt plus2pt}%%
    \list{}{%
%    \listparindent18pt\raggedright
        \leftmargin20pt
        \itemindent-10pt
    \parsep0pt
%    \itemsep3pt
    \partopsep0pt}
    \def\makelable##1{##1}%
}{%
    \endlist\par\addvspace{12pt plus2pt}%
}%
%
\def\webadd#1{\par\addvspace{3pt}{\fontfamily{Courier}\fontsize{9}{12}\selectfont\hskip18pt#1}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% itemize (bullet) %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% define \itemizeiTopSep if it is different from \topsepi
%
\newdimen\itemleftmargin \itemleftmargin1pt%\z@ % won't be active in enumerate
\newdimen\itemlabelsep \itemlabelsep5.5\p@ % labelsep in itemize for all levels
%
\def\itemize{%
  \ifnum \@itemdepth >3 \@toodeep\else\ifnum\@enumdepth>0\@itemdepth\@enumdepth\fi \advance\@itemdepth \@ne
     \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
     \list{\csname\@itemitem\endcsname}%
        {%
        \topsep\csname abovelistskip\romannumeral\the\@itemdepth\endcsname%\fi
        \itemsep\csname itemsep\romannumeral\the\@itemdepth\endcsname
        \labelsep5pt%\itemlabelsep
        \listfont
        \setbox\tempbox\hbox{\csname\@itemitem\endcsname}
        \csname leftmargin\romannumeral\the\@itemdepth\endcsname\wd\tempbox
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemlabelsep
        \ifnum\@itemdepth=\@ne\ifnum\@enumdepth=0 % if not inside an enumerate
           \rightmargin\listtextrightmargin
           \advance\rightmargin\rightskip
           \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmargini\leftskip
        \else
        \advance\csname leftmargin\romannumeral\the\@itemdepth\endcsname\itemleftmargin
           \advance\leftmarginii\leftskip
        \fi\fi
        \leftmargin\csname leftmargin\romannumeral\the\@itemdepth\endcsname
        \def\makelabel##1{\hss\llap{##1}}}
   \fi}%
%
\def\enditemize{%
\@topsepadd\csname belowlistskip\romannumeral\the\@itemdepth\endcsname%\fi
\endlist}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% enumroman (i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newenvironment{bulletlist}{%
   \abovelistskipii9pt plus2pt%
   \belowlistskipii9pt plus2pt%
   \renewcommand\labelitemi{\textbullet}%
   \renewcommand\labelitemii{\textbullet}%
   \renewcommand\labelitemiii{\textbullet}%
   \renewcommand\labelitemiv{\textbullet}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{hyphenlist}{%
   \renewcommand\labelitemi{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemii{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemiii{\normalfont\bfseries \textendash}%
   \renewcommand\labelitemiv{\normalfont\bfseries \textendash}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{asterisklist}{%
   \renewcommand\labelitemi{\textasteriskcentered}%
   \renewcommand\labelitemii{\textasteriskcentered}%
   \renewcommand\labelitemiii{\textasteriskcentered}%
   \renewcommand\labelitemiv{\textasteriskcentered}%
   \begin{itemize}
}{%
   \end{itemize}}
%
\newenvironment{caselist}{%
   \def\textlistlabel{Case\ }%
   \begin{enumerate}
}{%
   \end{enumerate}
}
%
\newenvironment{steplist}{%
   \def\textlistlabel{Step\ }%
   \begin{enumerate}
}{%
   \end{enumerate}
}
%
\newenvironment{labellist}{%
   \begin{description}
}{%
   \end{description}
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% npara %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{ncount}
\newenvironment{npara}{\par%
   \normalfont
   \setcounter{ncount}{0}%
   \def\item{\par\addtocounter{ncount}{1}\arabic{ncount}.\hskip2.6pt\ignorespaces}%
}{%
   \par}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Macro 2: TheoremLike Env. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
\newcounter{theorem}
\renewcommand\thetheorem{\thechapter.\arabic{theorem}}
%
\def\newtheorem#1{\@starthm{#1}\@ifnextchar[{\@othm{#1}}{\@nthm{#1}}}
\def\@starthm#1{\@namedef{#1*}{\@beginstarthm{#1}}%
	\@namedef{end#1*}{\@endstarthm{#1}}}
%
\def\@starthm#1{\@namedef{#1*}{\@beginstarthm{#1}}%
	\@namedef{end#1*}{\@endstarthm}}
%
\def\@beginstarthm#1{\@ifnextchar[{\@opargbeginstarthm{#1}}{\@opargbeginstarthm{#1}[]}}
%
\def\@opargbeginstarthm#1[#2]{\list{}{\topsep12pt plus2pt
 \def\@tempaa{#2}%
	\labelwidth0pt\labelsep9.5pt
	\leftmargin0pt
	\listparindent\parindent
 \edef\@Tempa{\csname #1name\endcsname}
	\def\makelabel##1{##1.}}
      \item[\hskip \labelsep{\csname#1headfont\endcsname \@Tempa{}\ifx\@tempaa\@empty\else\ (#2)\fi}]\csname#1font\endcsname}
%
\def\@endstarthm{\endlist}
%
\def\@xnthm#1#2[#3]{%
  \expandafter\@ifdefinable\csname #1\endcsname
    {\@definecounter{#1}\@newctr{#1}[#3]%
     \expandafter\xdef\csname the#1\endcsname{%
       \expandafter\noexpand\csname the#3\endcsname \@thmcountersep
          \@thmcounter{#1}}%
     \global\@namedef{#1}{\@thm{#1}{#2}{#1}}%
     \global\@namedef{end#1}{\@endtheorem}}}
\def\@ynthm#1#2{%
  \expandafter\@ifdefinable\csname #1\endcsname
    {\@definecounter{#1}%
     \expandafter\xdef\csname the#1\endcsname{\@thmcounter{#1}}%
     \global\@namedef{#1}{\@thm{#1}{#2}{#1}}%
     \global\@namedef{end#1}{\@endtheorem}}}
\def\@othm#1[#2]#3{%
  \@ifundefined{c@#2}{\@nocounterr{#2}}%
    {\expandafter\@ifdefinable\csname #1\endcsname
    {\global\@namedef{the#1}{\@nameuse{the#2}}%
  \global\@namedef{#1}{\@thm{#2}{#3}{#1}}%
  \global\@namedef{end#1}{\@endtheorem}}}}
\def\@thm#1#2#3{%
  \refstepcounter{#1}%
  \@ifnextchar[{\@ythm{#1}{#2}{#3}}{\@xthm{#1}{#2}{#3}}}
%
\def\@xthm#1#2#3{%
  \@opargbegintheorem{#2}{\csname the#1\endcsname}{}{#3}\ignorespaces}
\def\@ythm#1#2#3[#4]{%
  \@opargbegintheorem{#2}{\csname the#1\endcsname}{#4}{#3}\ignorespaces}
%
\def\@endtheorem{\par\addvspace{12pt plus2pt}\endlist}%
%
\def\@opargbegintheorem#1#2#3#4{\par\addvspace{12pt plus2pt}%
   \def\@tempa{#3}%
   \noindent{\csname #4headfont\endcsname\ifx\@tempa\empty#1\ #2.\hskip9.7pt\else#1\ #2\ {\reset@font\itshape#3}.\hskip9.7pt\fi}%
   \csname #4font\endcsname\ignorespaces}%
%
\def\@endtheorem{\par\addvspace{12pt plus2pt}}%
%
\def\examplename{Example}
\newtheorem{example}{\examplename}[chapter]
\def\theexample{\arabic{example}}
\def\examplefont{}
\def\exampleheadfont{\itshape}
%
\def\theoremname{Theorem}
\newtheorem{theorem}{\theoremname}[section]
\def\thetheorem{\arabic{theorem}}
\def\theoremfont{\itshape}
\def\theoremheadfont{\bfseries}
%
\def\lemmaname{Lemma}
\newtheorem{lemma}{\lemmaname}[section]
\def\thelemma{\arabic{lemma}}
\def\lemmafont{\itshape}
\def\lemmaheadfont{\bfseries}
%
\def\propositionname{Proposition}
\newtheorem{proposition}{\propositionname}[section]
\def\theproposition{\arabic{proposition}}
\def\propositionfont{\itshape}
\def\propositionheadfont{\bfseries}
%
\def\corollaryname{Corollary}
\newtheorem{corollary}{\corollaryname}[section]
\def\thecorollary{\arabic{corollary}}
\def\corollaryfont{\itshape}
\def\corollaryheadfont{\bfseries}
%
\def\definitionname{Def{}inition}
\newtheorem{definition}{\definitionname}[section]
\def\thedefinition{\arabic{definition}}
\def\definitionfont{}
\def\definitionheadfont{\bfseries}
%
\def\notationname{Notation}
\newtheorem{notation}{\notationname}[section]
\def\thenotation{\arabic{notation}}
\def\notationfont{\itshape}
\def\notationheadfont{\bfseries}
%
\def\assumptionname{Assumption}
\newtheorem{assumption}{\assumptionname}[section]
\def\theassumption{\arabic{assumption}}
\def\assumptionfont{}
\def\assumptionheadfont{\bfseries}
%
\def\remarkname{Remark}
\newtheorem{remark}{\remarkname}[section]
\def\theremark{\arabic{remark}}
\def\remarkfont{}
\def\remarkheadfont{\bfseries}
%
\newif\ifQEDlogo
\QEDlogofalse
\font\QEDlogofont=msam10 at 9pt
\def\QEDlogo{%
        \ifQEDlogo
 	\else
	  \ifvmode
	    \vskip-\lastskip\hspace*{\fill}%
          \else
  	    \unskip\hskip5pt\ignorespaces
          \fi
    	  \hfill{QED}%
          \global\QEDlogofalse
        \fi}
%
%%%%%% Proof %%%%%%
%
\def\prooffont{\reset@font\fontfamily{Fruti}\fontsize{9.5}{12}\itshape\selectfont}
\font\QEDlogofont=msam10 at 10pt
\def\QEDlogo{\hbox{\QEDlogofont\char'003}}
\def\QEDblogo{\hbox{\QEDlogofont\char'004}}
\newif\ifnologo\nologofalse
\newif\iflogo
\newif\ifblogo\blogofalse
\def\prbox{\unskip\hskip4.5pt\ignorespaces{\fboxsep0pt\fbox{\vrule width5pt height5pt depth0pt}}\global\logofalse}
\newenvironment{proof}{\par%
\leftskip18pt\rightskip\leftskip
			\par%
   \global\QEDlogofalse
   \@ifnextchar[{\@proof}{\@proof[{\prooffont Proof}]}%
   }{%
   \ifblogo\hskip1.2pt
   \QEDblogo
   \else
   \ifnologo
   \else
   \fi\fi
   \par%
   }
%
\def\@proof[#1]{\noindent{\prooffont#1\hskip4.8pt}\par\noindent\ignorespaces}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End TheoremLike Env. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\marginfont{\fontsize{9}{10}\selectfont\bfseries\raggedright}
\newcommand{\comment}[1]{%
 \setlength\overfullrule{0pt}%
 \ifodd\c@page
    \marginpar{%
    \parbox{8pc}{\par\vskip-1.5pt\marginfont please check.\par}}%
    \fboxsep.1pt\fbox{#1}%
 \else
    \marginpar{%
    \parbox{8pc}{\par\vskip-1.5pt\marginfont please check.\par}}
    \fboxsep.1pt\fbox{#1}\fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{topnumber}{5}
\renewcommand\topfraction{.9}
\setcounter{bottomnumber}{5}
\renewcommand\bottomfraction{.9}
\setcounter{totalnumber}{10}
\renewcommand\textfraction{.09}
\renewcommand\floatpagefraction{.901}
\setcounter{dbltopnumber}{1}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.901}
%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{9.5\p@}
\setlength\belowcaptionskip{.5\p@}
%
\def\FigName{figure}
\def\BoxName{boxenv}
%
\long\def\@makecaption#1#2{%
    \ifx\FigName\@captype
      \vskip\abovecaptionskip
      \@makefigurecaption{#1}{#2}%
    \else
      \ifx\BoxName\@captype
        \@makeboxenvcaption{#1}{#2}%
        \vskip\belowcaptionskip
      \else%
        \@maketablecaption{#1}{#2}%
        \vskip\belowcaptionskip%
      \fi%
    \fi%
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{figure}[chapter]
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{bp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename\ \thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
%
%
\def\@makefigurecaption#1#2{\figcaptionfont%
    \setbox\tempbox\hbox{\figcaptionfont{\figcaptionnumfont #1.}\enspace #2\par}%
    \ifdim\figwidth<19.6pc% 
        \ifdim\wd\tempbox>\columnwidth%
          {\figcaptionfont{\figcaptionnumfont#1.}\enspace #2\vphantom{pqg}\par}%
        \else
         \centerline{\figcaptionfont{\figcaptionnumfont#1.}\enspace #2}%
        \fi%
    \else%
        \ifdim\wd\tempbox>\textwidth%
          {\figcaptionfont{\figcaptionnumfont#1.}\enspace #2\vphantom{pqg}\par}%
        \else
         \centerline{\figcaptionfont{\figcaptionnumfont#1.}\enspace #2}%
        \fi%
  \fi%
}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{table}[chapter]
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{bp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\ \thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
%
\newenvironment{Table}
               {\par\addvspace{12pt plus2pt}\def\@captype{table}}
               {\par\addvspace{12pt plus2pt}}
%
\newenvironment{Table*}
               {\par\addvspace{12pt plus2pt}\def\@captype{table}}
               {\par\addvspace{12pt plus2pt}}
%
\def\@maketablecaption#1#2{%
    {\tablecaptionfont\noindent{\tablecaptionnumfont#1.}\enspace #2\vphantom{y}\par}}%
%
\def\tbl#1#2{{\def\arraystretch{1.5}%%
  \setbox\tempbox\hbox{\tablefont#2}%
  \tabledim\hsize
  \advance\tabledim-\wd\tempbox
  \ifdim\tabledim>0pt
    \divide\tabledim2
  \else
    \global\tabledim0pt
  \fi
  \global\tableleftskip\tabledim plus1fill%
  \global\tablerightskip\tabledim plus1fill%
  \global\tabnoteleftskip\tabledim% plus1fill%
  \global\tabnoterightskip\tabledim plus1fill%
  
  \caption{#1}%
  {\centerline{\box\tempbox}}}}%
%
\def\x{@{\extracolsep{\fill}}}
\def\toprule{\\[-10pt]\Hline\\[-12pt]\Hline\\[-9pt]}%
\def\colrule{\\[-12pt]\Hline\\[-9pt]}%
\def\botrule{%%
\\[-11pt]\Hline}%
\def\crule#1{\\[-5pt]#1\\[-4.5pt]}
%
\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height .25pt \futurelet%\@height \arrayrulewidth
   \@tempa\@xhline}
%
\newenvironment{tabnote}{\par\addvspace{5pt}%
  \noindent\tabnotefont\ignorespaces%
  }{\par}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Box         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{boxenv}[chapter]
\renewcommand\theboxenv{\@arabic\c@boxenv}
\def\fps@boxenv{tbp}
\def\ftype@boxenv{2}
\def\ext@boxenv{lob}
\newcommand\boxenvname{Box}
\def\fnum@boxenv{\boxenvname\ \theboxenv}
%
\newenvironment{boxenv}
               {\@float{boxenv}}
               {\end@float}
%
\def\@makeboxenvcaption#1#2{%
    {\boxcaptionfont\noindent{\boxcaptionnumfont#1.}\hskip3.5pt#2\vphantom{y}\par}}%
%
\def\boxdata#1#2{%
  \caption{#1}%
  {\boxfont\FrameSep6pt\FrameRule.5pt\topsep0pt\begin{framed}#2\end{framed}}%
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%\def\arraystretch{2}%
\setlength\arraycolsep{1.5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
%
\setlength\columnsep{18\p@}
\setlength\columnseprule{0\p@}
%
\def\tstrut{\vrule width0pt depth4.65pt}
\def\bstrut{\vrule width0pt height9pt}
\jot=3pt
%
\def\text#1{\mathchoice
  {\hbox{\fontsize{\tf@size}{\tf@size}\selectfont#1}}%
  {\hbox{\fontsize{\tf@size}{\tf@size}\selectfont#1}}%
  {\hbox{\fontsize{\sf@size}{\sf@size}\selectfont#1}}%
  {\hbox{\fontsize{\ssf@size}{\ssf@size}\selectfont#1}}}
%
\mathindent1em%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Math %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Footnote %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule width6.5pc height.25pt%
  \kern 5\p@}
%
\newdimen\@footmax
\def\footmax#1{%
  \setbox\tempbox\hbox{\footnotesize$^{#1}$\hskip.5em}%
  \global\@footmax\wd\tempbox}
%
\footmax{1}
%
\long\def\@makefntext#1{\parindent9pt%
    {\leavevmode\let\doPar\relax{\@makefnmark}#1}%
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Footnote %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Page styles  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
\def\jvol#1{\gdef\@jvol{#1}}
\def\jnum#1{\gdef\@jnum{#1}}
\def\jmonth#1{\gdef\@jmonth{#1}}
\def\jyear#1{\gdef\@jyear{#1}}
\def\jname#1{\gdef\@jname{#1}}
%
\def\@jvol{}
\def\@jnum{}
\def\@jmonth{}
\def\@jyear{}
\def\@jname{J. Audio Eng. Sco.}
%
\def\oddfolio{{\foliofont\thepage}}
\def\evenfolio{{\foliofont\thepage}}
%
\def\ps@headings{%
    \def\@evenhead{{\rhfont\leftmark}\hfil{\rhfont PAPERS}}%
    \def\@oddhead{{\rhfont PAPERS}\hfil{\rhfont\rightmark}}%
    \def\@evenfoot{\hbox to \textwidth{\rhfont\evenfolio\hfill \@jname, Vol.\ \@jvol,\ No.\ \@jnum,\ \@jyear\ \@jmonth}}%
    \def\@oddfoot{\hbox to \textwidth{\rhfont \@jname, Vol.\ \@jvol,\ No.\ \@jnum,\ \@jyear\ \@jmonth\hfill\oddfolio}}%
    \let\@mkboth\markboth
    \def\chaptermark##1{\markboth{##1}{##1}}%
}
%
\def\ps@empty{%
      \def\@oddfoot{}%
      \let\@evenfoot\@oddfoot
      \def\@evenhead{}%
      \def\@oddhead{}%
      \let\@mkboth\markboth
      \let\chaptermark\@gobble
      \let\sectionmark\@gobble}
%
\def\Copyright{\raise.2ex\hbox{\fontsize{9}{9}\selectfont{\ooalign{\hfil\raise-.05ex\hbox{\kern-.3pt\fontsize{7}{8}\selectfont C}%
                \hfil\crcr \mathhexbox20D}}}}
%
\def\bottomfoot#1{\gdef\bottom@foot{#1}}
%
\def\ps@plain{\let\@mkboth\@gobbletwo%
       \def\@evenhead{\hbox to \textwidth{\hfill\rhfont PAPERS}}%
       \def\@oddhead{\hbox to \textwidth{\rhfont PAPERS\hfill}}%
       \def\@evenfoot{\hbox to \textwidth{\rhfont\thepage\hfill \@jname, Vol.\ \@jvol,\ No.\ \@jnum,\ \@jyear\ \@jmonth}}%
       \def\@oddfoot{\hbox to \textwidth{\rhfont \@jname, Vol.\ \@jvol,\ No.\ \@jnum,\ \@jyear\ \@jmonth\hfill\thepage}}%
}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Page styles  %%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Bibliography Section   %%%%%%%%%%%%%%%%%%
%
\newdimen\bibitemsep
\newdimen\bibparsep
\newdimen\biblabelsep
%
%%%%%%%% Unnumberd References %%%%%%%%%%
%
\setlength\bibitemsep  {0pt}
\setlength\bibparsep  {0pt}
%%%%%%%% Numbered Within Squre Bracket %%%%%%%%%%%
\newif\ifnumsqure
\global\numsqurefalse
\setlength\biblabelsep  {4.7pt}
%%%%%%%% Numbered With Dot %%%%%%%%%%%%%%%
%
\newenvironment{thebibliography}[1]{\overfullrule0pt%
     \bibliofont%     
     \global\labelsepi0pt
     \global\LabelSep0pt
     \global\leftmargini0pt
     \def\@tempa{#1}%
     \section{\bibname}
     \ifx\@tempa\@empty
         \list{}%
           {\labelwidth0pt\labelsep0pt
           \leftmargin0pt%
           \itemsep\bibitemsep
           \parsep\bibparsep
           \usecounter{enumiv}%
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\arabic{enumiv}}}%
      \else
          \setbox\tempbox\hbox{[\@tempa]}% \@tempa.
          \tempdimen\wd\tempbox
          \ifnumsqure
          \else
          \def\@biblabel##1{{[##1]}}\fi%
           \list{\@biblabel{\arabic{enumiv}}}%
           {%
	     \leftskip0pt
           \leftmargin.001pt
           \labelsep\biblabelsep
	   \listparindent0pt
           \itemindent4.5pt
           \advance\itemindent\biblabelsep
           \itemsep\bibitemsep
           \parsep\bibparsep
           \usecounter{enumiv}%
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\arabic{enumiv}}}%
     \fi
     \sfcode`\.=\@m}
     {\def\@noitemerr
     {\@latex@warning{Empty `thebibliography' environment}}%
     \endlist%
     \reset@font}
%
\newcommand\newblock{}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End Bibliography Section %%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg {2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{2}
%
\newcommand\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
       \nobreak
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}
    \endgroup
  \fi}
\newcommand\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}
\newcommand\l@section      {\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\l@subsection   {\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand\l@paragraph    {\@dottedtocline{4}{10em}{5em}}
\newcommand\l@subparagraph {\@dottedtocline{5}{12em}{6em}}
\let\l@table\l@figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End FM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Appendix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand\appendix{\par
  \setcounter{chapter}{0}%
\renewcommand\thesection{\Alph{section}}
  \setcounter{section}{1}%
  \setcounter{equation}{0}%
  \let\section\subsection
  \renewcommand\thechapter{\Alph{chapter}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Index  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\alphabet#1{{\bf#1}\par\vskip.5\baselineskip}
%
\newenvironment{theindex}{%
      \columnseprule \z@
      \columnsep 12\p@
      \indexfont
      \begin{multicols}{4}
      \parindent\z@
}{\end{multicols}%
      \par}
\def\inditem#1#2{{\indexfont #1\hspace*{5pt}#2\par}}
\newcommand\@idxitem  {\par\hangindent3em}
\newcommand\subitem   {\par\hangindent3em\hspace*{1em}}
\newcommand\subsubitem{\par\hangindent3em\hspace*{2em}}
\newcommand\indexspace{\par \vskip\baselineskip}
%
%%%%%%%%%%%%%%%%%%%% File Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\document{\endgroup
  \ifx\@unusedoptionlist\@empty\else
    \@latex@warning@no@line{Unused global option(s):^^J%
            \@spaces[\@unusedoptionlist]}%
  \fi
  \@colht\textheight
  \@colroom\textheight \vsize\textheight
  \columnwidth\textwidth
  \@clubpenalty\clubpenalty
  \if@twocolumn
    \advance\columnwidth -\columnsep
    \divide\columnwidth\tw@ \hsize\columnwidth \@firstcolumntrue
  \fi
  \hsize\columnwidth \linewidth\hsize
  \begingroup\@floatplacement\@dblfloatplacement
    \makeatletter\let\@writefile\@gobbletwo
    \global \let \@multiplelabels \relax
    \@input{\jobname.aux}%
  \endgroup
  \if@filesw
    \immediate\openout\@mainaux\jobname.aux
    \immediate\write\@mainaux{\relax}%
  \fi
  \process@table
  \let\glb@currsize\@empty  %% Force math initialisation.
  \normalsize
  \everypar{}%
  \@noskipsecfalse
  \let \@refundefined \relax
  \let\AtBeginDocument\@firstofone
  \@begindocumenthook
  \ifdim\topskip<1sp\global\topskip 1sp\relax\fi
  \global\@maxdepth\maxdepth
  \global\let\@begindocumenthook\@undefined
  \ifx\@listfiles\@undefined
    \global\let\@filelist\relax
    \global\let\@addtofilelist\@gobble
  \fi
  \gdef\do##1{\global\let ##1\@notprerr}%
  \@preamblecmds
  \global\let \@nodocument \relax
  \global\let\do\noexpand
  \ignorespaces}
%
\def\cropmarks{\nointerlineskip\vbox to 0pt{\offinterlineskip\vss
  \hbox to \trimwidth{\vbox to 0pt{ }\hfill}}}%

\def\croprules{%
		\nointerlineskip\vbox to 0pt{\offinterlineskip
  \vskip4pt\cropmarks
		\vskip\topmargin
  \hbox to \trimwidth{\hskip\@themargin\hbox to \textwidth{\hss%
      \vrule height0pt depth\typeheight width.1pt%
      \vrule height.1pt depth0pt width\textwidth
      \vrule height0pt depth\typeheight width.1pt\hss}}%
  \hbox to \trimwidth{\hskip\@themargin\hbox to \textwidth{\hss%
      \vrule height.1pt depth0pt width\textwidth\hss}\hss}
}}
%
%%%%%%%%%%%%%%%%%%% End output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Routine  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%
\def\@outputpage{%
\begingroup           % the \endgroup is put in by \aftergroup
  \let \protect \noexpand
  \@resetactivechars
  \@parboxrestore
  \shipout \vbox{%
    \set@typeset@protect
    \aftergroup \endgroup
    \aftergroup \set@typeset@protect
                                % correct? or just restore by ending
                                % the group?
  \if@specialpage
    \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
  \fi
  \if@twoside
    \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
         \let\@themargin\oddsidemargin
    \else \let\@thehead\@evenhead
       \let\@thefoot\@evenfoot \let\@themargin\evensidemargin
    \fi
  \fi
  \reset@font
  \normalsize
  \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
    \@begindvi
    \vskip \topmargin
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@hbox
          \normalcolor
          \hb@xt@\textwidth {%
            \let \label \@gobble
            \let \index \@gobble
            \let \glossary \@gobble %% 21 Jun 91
            \@thehead
            }%
        \color@endbox
        }%                        %% 22 Feb 87
      \dp\@tempboxa \z@
      \box\@tempboxa
      \vskip \headsep
      \box\@outputbox
      \baselineskip \footskip
      \color@hbox
        \normalcolor
        \hb@xt@\textwidth{%
          \let \label \@gobble
          \let \index \@gobble      %% 22 Feb 87
          \let \glossary \@gobble   %% 21 Jun 91
          \@thefoot
          }%
      \color@endbox
      }%
    }%
\global \@colht \textheight
\stepcounter{page}%
\let\firstmark\botmark
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Output Routine  %%%%%%%%%%%%%%%%%%%%%%
%
%============================================
\newcount\authorcount
\newcount\affcount
\newcount\tempcount
%
\def\thanks#1{\footnotemark
    \protected@xdef\@thanks{\@thanks
        \protect\footnotetext[\the\c@footnote]{#1}}%
}
%
\def\equation{\@beginparpenalty\predisplaypenalty
  \@endparpenalty\postdisplaypenalty
		\refstepcounter{equation}\trivlist \item[]\leavevmode
  \hbox to\linewidth\bgroup $\m@th% $ TO MAKE DOLLAR NESTING OK
  \displaystyle
\mathindent10pt\hskip\mathindent}

\def\endequation{$\hfil
           \displaywidth\linewidth\@eqnnum\egroup \endtrivlist}
% \newenvironment{equation*}{\relax
%                 \ifmmode\@badmath
%                 \else
%                   \begin{trivlist}%
%                     \@beginparpenalty\predisplaypenalty
%                     \@endparpenalty\postdisplaypenalty
%                     \item[]\leavevmode
%                     \hb@xt@\linewidth\bgroup $\m@th\displaystyle %$
%                       \mathindent10pt\hskip\mathindent\bgroup
%                 \fi}{\relax
%                 \ifmmode
%                       \egroup $%\hfil% $
%                     \egroup
%                   \end{trivlist}%
%                 \else \@badmath
%                 \fi}
%
%
\def\eqnarray{\stepcounter{equation}%
\def\@currentlabel{\p@equation\theequation}%
\global\@eqnswtrue\m@th%
\mathindent10pt%
\global\@eqcnt\z@\tabskip\mathindent%
\let\\\@eqncr%
$$\everycr{}\halign to\displaywidth\bgroup\hskip\@centering%
$\displaystyle\tabskip\z@skip{##}$\@eqnsel&\global\@eqcnt\@ne%
\hskip\tw@\arraycolsep\hfil${##}$\hfil%
&\global\@eqcnt\tw@\hskip\tw@\arraycolsep%
$\displaystyle{##}$\hfil%
\tabskip\@centering&\global\@eqcnt\thr@@%
\hbox to\z@\bgroup\hss##\egroup\tabskip\z@skip\cr}
%
% \DeclareRobustCommand{\[}{\begin{equation*}}
% \DeclareRobustCommand{\]}{\end{equation*}}
%%
\def\onecolumn{%
  %\clearpage
  \global\columnwidth\textwidth
  \global\hsize\columnwidth
  \global\linewidth\columnwidth
  \global\@twocolumnfalse
  \col@number \@ne
  \@floatplacement}

%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%% SL 001202 - begin
%
\def\@outputdblcol{%
  \if@firstcolumn
    \global \@firstcolumnfalse
    \global \setbox\@leftcolumn \box\@outputbox
    \ifswitchcols\global\switchcolsfalse
       \ifdim\@partialpageht>0pt\advance\@partialpageht-.5pt\fi
       \vglue\@partialpageht
    \fi
  \else
    \global \@firstcolumntrue
    \setbox\@outputbox \vbox {%
                         \hbox to\textwidth {%
                           \hbox to\columnwidth {%
                             \box\@leftcolumn \hss}%
                           \hfil
                           \vrule \@width\columnseprule
                           \hfil
                           \hbox to\columnwidth {%
                             \box\@outputbox \hss}%
                                             }%
                              }%
    \@combinedblfloats
    \@outputpage
    \begingroup
      \@dblfloatplacement
      \@startdblcolumn
      \@whilesw\if@fcolmade \fi
        {\@outputpage
         \@startdblcolumn}%
    \endgroup
  \fi
}
%
\let\save@clearpage\clearpage
\let\save@outputpage\@outputpage
\newdimen\@partialpageht\@partialpageht10pc
\newif\ifswitchcols
\newbox\store@outputbox
\newdimen\storedboxheight
%
\def\store@outputpage{%
  \global\setbox\store@outputbox \vbox to \storedboxheight{%
      \leftline{\box\@outputbox}\vss}%
   \noindent%\fbox
   {\box\store@outputbox}%
   \par%\vrule width\textwidth height1pt
   \global\let\@outputpage\save@outputpage}
%
\def\balancepage{\@ifnextchar[{\balancecolpagebreak}{\balancecolpagebreak[0pt]}}
\def\balancecolpagebreak[#1]{\ifvmode\vadjust\fi{\vskip-\baselineskip
                             \leftline{\vrule width\typewidth height\draftrule depth0pt}%
                             \expandafter\ifdim#1=0pt\vfill\else\vskip#1\fi}\pagebreak}
%
\def\balanceandcontinue{%
   \par
   \ifdim\pagegoal=\maxdimen
   \else
      \global\let\@outputpage\store@outputpage
      \ifdim\@colht=\pagegoal
         \global\storedboxheight\pagetotal
      \else
         \global\storedboxheight\@colht
         \global\advance\storedboxheight-\pagegoal
         \global\advance\storedboxheight\pagetotal
      \fi
   \fi
   \onecolumn
   \global\switchcolstrue}
%
\def\TwoColumn{\@ifnextchar[{\@TwoColumn}{\@TwoColumn[]}}
\def\@TwoColumn[#1]{%
   \balanceandcontinue
   \let\clearpage\relax
   \def\@tempa{#1}%
   \ifx\@tempa\@empty
      \twocolumn
   \else
      #1\global\@partialpageht\pagetotal\twocolumn
   \fi
   \let\clearpage\save@clearpage
}
%
\font\aulcir = lcircle10 at 20pt
\newcommand\aubulls{\raise1.5pt\hbox{\aulcir\char'162}}
\def\autextbullet{\leavevmode\raise1pt\hbox{\aubulls}}
%
\newdimen\nommaxdim
\newbox\nombox
\def\nomentry#1#2{{\par\@hangfrom{\hbox to \nommaxdim{\hfill{}#1} $=$ }#2\par}}
%
\newenvironment{nomenclature}[1][\relax]{\par%
                               \ifx#1\relax\else\setbox\nombox\hbox{#1}\global\nommaxdim\wd\nombox\fi}%
                               {\par}
%
\RequirePackage{cuted}%
\stripsep0pt
%%%%%%%%Modified cuted macro
\gdef\@viperoutput{%
               \global \setbox\@holdpg \vbox{%
		\vrule height5pt width0pt depth0pt%
                      \unvbox\@holdpg
                      \unvbox\@cclv
                      \setbox\@tempboxa \lastbox
                      \unskip
                                     }%
               \@pagedp \dp\@holdpg
               \@pageht \ht\@holdpg
               \ifvoid\@viper\else
                   \advance \@pageht \@pagedp
                   \ifvoid\footins \else
                     \advance \@pageht \ht\footins
                     \advance \@pageht \skip\footins
                     \advance \@pageht \dp\footins
                   \fi
                   \ifvbox \@kludgeins
                     \ifdim \wd\@kludgeins=\z@
                       \advance \@pageht \ht\@kludgeins
                     \fi
                   \fi
                   \@reinserts
                   \@addviper
               \fi
}
%%%%%%%%%
\newcount\calccount
\newcount\biocount
\newcount\printbiocount
\def\biography#1#2#3{%
\global\advance\biocount\@ne
\expandafter\gdef\csname bioname\romannumeral\biocount\endcsname{#1}%
\expandafter\gdef\csname biofig\romannumeral\biocount\endcsname{#2}%
\expandafter\gdef\csname bioaff\romannumeral\biocount\endcsname{#3}}

\def\printbiography{%%
\authorcount\biocount%
\ifnum\biocount=\tw@ \tabcolsep5pc\fi%
\affcount=\z@%
\printbiocount\biocount%
\tempcount=\z@%
\begin{strip}\par\addvspace{\baselineskip}
\noindent\vrule height1pt width\textwidth\par\addvspace{\baselineskip}
{\sectionfont\centering THE \ifnum\biocount=\@ne AUTHOR\else AUTHORS\fi\par\addvspace{\baselineskip}}%
\begingroup%
\loop%
\ifnum\printbiocount>0%
  \calccount\printbiocount%
  \advance\calccount-4
  \global\advance\printbiocount-4%
      \centering\@whilenum\biocount>\calccount%
    \do{%
       \ifnum\biocount>0%
       \begin{tabular}{c}%
       \global\advance\tempcount\@ne%
       \edef\temp{\csname biofig\romannumeral\tempcount\endcsname}%
       \IfFileExists{\temp}{\includegraphics[width=8pc]{\temp}}{\vrule width8pc height10pc}\\%
       \csname bioname\romannumeral\tempcount\endcsname%
       \end{tabular}%
        \fi%
       \global\advance\biocount\m@ne%
    }%
  \\[1.5\baselineskip]%
\repeat%
\endgroup%
\end{strip}
\@whilenum\authorcount>\z@%%
\do{%%
  \advance\affcount\@ne%%
  \advance\authorcount\m@ne%%
  \ifnum\affcount>\@ne \centerline{\autextbullet}\vskip6pt\fi%%
  \par\csname bioaff\romannumeral\affcount\endcsname\par\addvspace{.5\baselineskip}%%
}%
\begin{strip}\end{strip}
\noindent\vrule height1pt width\textwidth\par\addvspace{\baselineskip}
}%
\AtEndDocument{\ifauthorgd\else\printbiography\fi}%

%%%%%%%%%%%%%%%%%%%%%%%%% Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\RequirePackage{graphicx}
%
\medmuskip=4mu 
\thickmuskip=5mu 
%
\pagestyle{headings}
\pagenumbering{arabic}

\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
%
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty

%%
\tolerance 1414
\hbadness 1414
\emergencystretch 1.5em
\hfuzz 0.3pt
\widowpenalty=10000
\clubpenalty10000
\vfuzz \hfuzz
\frenchspacing
%
\voffset-6pc
\hoffset-6pc


\endinput
%
