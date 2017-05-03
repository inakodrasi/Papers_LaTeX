% This is aes60i.cls for the AES 60th International Conference
% A skeleton class to load the main class aes.cls
\ProvidesClass{aes60i}
         [2002/02/04 v0.01 AES 60th International Conference paper class]
% 60th International Conference, Leuven, Belgium, 2016 February 3--5
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{aes}}
\ProcessOptions
\def\AEScontype{Conference}%
\def\AESconname{60th International \AEScontype}%
\def\AEScondate{2016 February 3--5}%
\def\AESconplace{Leuven, Belgium}%
\def\AESconRHname{AES 60${}^{\mbox{\scriptsize TH}}$ INTERNATIONAL CONFERENCE}%
\def\AESconRHdate{2016 February 3--5}%
\def\AESconRHplace{Leuven, Belgium}%
\LoadClass[mathptmx]{aes}[2006/04/02]
\endinput

