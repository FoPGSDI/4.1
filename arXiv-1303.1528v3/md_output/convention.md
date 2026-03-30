# Notation & Conventions

## Custom Commands

| Command | Args | Expansion | Description |
|---------|------|-----------|-------------|
| `\note` | 1 | `[$\blacktriangleright$~\textbf{#1}~$\blacktriangleleft$]` |  |
| `\ch` | 0 | `[$\blacktriangleright$~{\bf{Checked!}}~$\blacktriangleleft$]~` |  |
| `\be` | 0 | `\begin{equation}` |  |
| `\ba` | 0 | `\begin{eqnarray}` |  |
| `\ee` | 0 | `\end{equation}` |  |
| `\ea` | 0 | `\end{eqnarray}` |  |
| `\bw` | 0 | `\begin{widetext}` |  |
| `\ew` | 0 | `\end{widetext}` |  |
| `\mb` | 1 | `\mbox{\boldmath $#1$}` |  |
| `\comment` | 1 | `{\bf [#1]}` |  |
| `\nn` | 0 | `\nonumber` |  |
| `\Schw` | 0 | `{\mbox{\tiny Schw}}` |  |
| `\hGR` | 0 | `h` |  |
| `\hDef` | 0 | `\mathfrak{h}` |  |
| `\ppE` | 0 | `{\mbox{\tiny ppE}}` |  |
| `\Kerr` | 0 | `{\mbox{\tiny Kerr}}` |  |
| `\ISCO` | 0 | `{\mbox{\tiny ISCO}}` |  |
| `\GR` | 0 | `{\mbox{\tiny GR}}` |  |
| `\ZM` | 0 | `{\mbox{\tiny ZM}}` |  |
| `\RW` | 0 | `{\mbox{\tiny RW}}` |  |
| `\IZ` | 0 | `{\mbox{\tiny IZ}}` |  |
| `\NZ` | 0 | `{\mbox{\tiny NZ}}` |  |
| `\FZ` | 0 | `{\mbox{\tiny FZ}}` |  |
| `\TT` | 0 | `{\mbox{\tiny TT}}` |  |
| `\MAT` | 0 | `{\mbox{\tiny mat}}` |  |
| `\STF` | 0 | `{\mbox{\tiny STF}}` |  |
| `\CPM` | 0 | `{\mbox{\tiny CPM}}` |  |
| `\CS` | 0 | `{\mbox{\tiny CS}}` |  |
| `\GW` | 0 | `{\mbox{\tiny GW}}` |  |
| `\AXIAL` | 0 | `{\mbox{\tiny Axial}}` |  |
| `\POLAR` | 0 | `{\mbox{\tiny Polar}}` |  |
| `\MIN` | 0 | `{\mbox{\tiny min}}` |  |
| `\INC` | 0 | `{\mbox{\tiny inc}}` |  |
| `\APO` | 0 | `{\mbox{\tiny apo}}` |  |
| `\PERI` | 0 | `{\mbox{\tiny peri}}` |  |
| `\pont` | 0 | `{\,^\ast\!}R\,R` |  |
| `\met` | 0 | `\mbox{g}` |  |
| `\metb` | 0 | `\mbox{\bf g}` |  |
| `\oone` | 0 | `{}^{\mbox{\tiny $(1)$}}` |  |
| `\otwo` | 0 | `{}^{\mbox{\tiny $(2)$}}` |  |
| `\lgw` | 0 | `\lambda^{}_{\GW}` |  |
| `\hgw` | 0 | `h^{}_{\GW}` |  |
| `\tamaE` | 1 | `\textcolor{blue}{\sout{#1}}` |  |
| `\change` | 1 | `\textcolor{red}{#1}` |  |
| `\kent` | 1 | `\textcolor{green}{\textbf{ #1}}` |  |
| `\pd` | 0 | `\partial` |  |
| `\cd` | 0 | `\nabla` |  |
| `\LevTen` | 0 | `\varepsilon` |  |
| `\JJ` | 4 | `J^{(#1)}_{#2,#3,#4}` |  |
| `\Jdef` | 4 | `\lim_{3\to #1}\pd^{(1)}_{#2}\pd^{(2)}_{#3}\pd^{(3)}_{#4}\mathcal{G}(ABC)` |  |
| `\plusonetotwo` | 0 | `+(1\leftrightarrow 2)` |  |
| `\plusitoj` | 0 | `+(i\leftrightarrow j)` |  |
| `\ttilde` | 1 | `\tilde{\tilde{#1}}` |  |
| `\Q` | 0 | `{\mbox{\tiny Q}}` |  |
| `\Def` | 0 | `{\mbox{\tiny Def}}` |  |
| `\boldh` | 0 | `\mathbf{h}` |  |
| `\cs` | 0 | `{\mbox{\tiny CS}}` |  |
| `\gr` | 0 | `{\mbox{\tiny GR}}` |  |
| `\INT` | 0 | `{\mbox{\tiny INT}}` |  |
| `\PP` | 0 | `{\mbox{\tiny PP}}` |  |
| `\BH` | 0 | `{\mbox{\tiny BH}}` |  |
| `\NPP` | 0 | `{\mbox{\tiny NPP}}` |  |
| `\self` | 0 | `{\mbox{\tiny self}}` |  |
| `\cross` | 0 | `{\mbox{\tiny cross}}` |  |
| `\NS` | 0 | `*` |  |
| `\ext` | 0 | `\mathrm{ext}` |  |
| `\inter` | 0 | `\mathrm{int}` |  |
| `\RNS` | 0 | `\mathcal{R}_*` |  |
| `\tid` | 0 | `\mathrm{(tid)}` |  |
| `\rot` | 0 | `\mathrm{(rot)}` |  |
| `\N` | 0 | `{\mbox{\tiny N}}` |  |
| `\mrm` | 0 | `\mathrm` |  |
| `\lb` | 0 | `\left(` |  |
| `\rb` | 0 | `\right)` |  |
| `\lcb` | 0 | `\left\{}
\newcommand{\rcb}{\right\}` |  |
| `\lsb` | 0 | `\left[` |  |
| `\rsb` | 0 | `\right]` |  |
| `\ld` | 0 | `\left.` |  |
| `\rd` | 0 | `\right.` |  |
| `\E` | 0 | `{\mbox{\tiny E}}` |  |
| `\accR` | 0 | `\mathscr{R}` |  |
| `\accS` | 0 | `\mathscr{S}` |  |
| `\accW` | 0 | `\mathscr{W}` |  |
| `\lambdabartid` | 0 | `\bar{\lambda}^\mathrm{(tid)}` |  |
| `\Ibar` | 0 | `\bar{I}` |  |
| `\Qbar` | 0 | `\bar{Q}` |  |
| `\lambdabarrot` | 0 | `\bar{\lambda}^\mathrm{(rot)}` |  |
| `\ave` | 1 | `\left< #1 \right>` |  |
| `\arraystretch` | 0 | `1.2` |  |
| `\arraystretch` | 0 | `1.2` |  |
