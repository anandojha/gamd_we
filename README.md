# gamd_we - Gaussian Accelerated Molecular Dynamics - Weighted Ensemble

## Installation and Setup Instructions :
* Make sure [anaconda3](https://www.anaconda.com/) is installed on the local machine. Go to the  [download](https://www.anaconda.com/products/individual) page of anaconda3 and install the latest version of anaconda3.
* Create a new conda environment with python = 3.6 and install the package with the following commands in the terminal: 
```bash
conda create -n gamdwe python=3.6
conda activate gamdwe
conda install -c conda-forge curl matplotlib openmm seaborn pandas 
conda install -c ambermd pytraj
conda install -c omnia ambertools

## GaMD :
\begin{itemize}
\item GaMD enhances the conformational sampling by adding a harmonic boost potential to smoothen the system potential energy surface.
\item Consider a system with N atoms at positions, $\vec{r}$ = \{$\vec{r\textsubscript{1}}$, $\vec{r\textsubscript{2}}$, $\vec{r\textsubscript{3}}$,  . . ,$\vec{r\textsubscript{N}}$\}
\item For a threshold energy, E and V($\vec{r}$) $<$ E, we add a boost potential, $\Delta$ V such that  V\textsuperscript{*}($\vec{r}$) = V($\vec{r}$) + $\Delta$ V
\item Two necessary and satisfying conditions: \par 
\begin{enumerate}
\item V\textsubscript{1}($\vec{r}$) < V\textsubscript{2}($\vec{r}$)  \implies {V\textsubscript{1}}\textsuperscript{*}(\vec{r}) < {V\textsubscript{2}}\textsuperscript{*}(\vec{r}) \implies E < \frac{1} {2} (V\textsubscript{1}(\vec{r}) + V\textsubscript{2}(\vec{r}))  + $\frac{1}$ / {k}
\item V\textsubscript{1}($\vec{r}$) < V\textsubscript{2}($\vec{r}$)  \implies {V\textsubscript{1}}\textsuperscript{*}(\vec{r}) -  {V\textsubscript{2}}\textsuperscript{*}(\vec{r}) < V\textsubscript{1}(\vec{r}) - V\textsubscript{2}(\vec{r})
\implies E > \frac{1} {2} (V\textsubscript{1}(\vec{r}) + V\textsubscript{2}(\vec{r})) \par 
\end{enumerate}

\item From above two relationships:  $\frac{1} {2}$ (V\textsubscript{1}($\vec{r}$) + V\textsubscript{2}($\vec{r}$)) $<$  E $<$ $\frac{1} {2}$ (V\textsubscript{1}($\vec{r}$) + V\textsubscript{2}($\vec{r}$))  +  $\frac{1}{k}$
\item  V\textsubscript{min} = Minimum PE and  V\textsubscript{max} = Maximum PE \par \implies V\textsubscript{min} $\leq$ V\textsubscript{1}($\vec{r}$) + V\textsubscript{2}($\vec{r}$) $\leq$ V\textsubscript{max} $\implies$ V\textsubscript{max} $\leq$  E $\leq$ V\textsubscript{min} + $\frac{1}{k}$ $\centernot\implies$ k $\leq$ $\frac {1} {V\textsubscript{max} - V\textsubscript{min}}$
\item k = k\textsubscript{0} ($\frac {1} {V\textsubscript{max} - V\textsubscript{min}}$) for 0 $<$ k\textsubscript{0} $\leq$  1

\item Case I: When E = V\textsubscript{max} (Lower Bound) \par 
k\textsubscript{0} = min (1.0, ($\frac{\sigma \textsubscript{0}}{\sigma \textsubscript{V}}$)($\frac{V\textsubscript{max} - V\textsubscript{min}}{V\textsubscript{max} - V\textsubscript{avg}}$)) \par 
where $\sigma\textsubscript{V}$ = standard deviation of potential energies,  and $\sigma\textsubscript{0}$ = user-defined upper limit 
\item Case II: When E = V\textsubscript{min} + $\frac{1}{k}$ (Upper Bound) \par 
k\textsubscript{0} $\geq$ (1 - $\frac{\sigma \textsubscript{0}}{\sigma \textsubscript{V}}$)($\frac{V\textsubscript{max} - V\textsubscript{min}}{V\textsubscript{avg} - V\textsubscript{min}}$) for 0 < k\textsubscript{0} $\leq$ 1 and E = V\textsubscript{max} for k\textsubscript{0} >1


\item p (A\textsubscript{j}) = $\frac{p\textsuperscript{*}(A\textsubscript{j}) < e \textsuperscript{ \beta\Delta V ($\vec${r})} >  }{\smashoperator{\sum_{j = 1}^{M}} < p\textsuperscript{*}(A\textsubscript{j}) < e \textsuperscript{ \beta\Delta V ($\vec${r})} >  _j}$
\par
where p\textsuperscript{*}(A\textsubscript{j}) = probability distribution along a certain coordinate, M = number of bins, \par  < e \textsuperscript{\beta\Delta V ($\vec${r})} >  = ensembled average boltzmann factor of $\Delta$ V for simulation frames in the j\textsubscript{th} bin, and $\beta$ = k\textsubscript{B} T


\item Cumulant Expansion: < e \textsuperscript{\beta\Delta V ($\vec${r})} >  =  e$^ {{\smashoperator{\sum_{k = 1}^{\infty}} \frac{\beta \textsuperscript{k}{}}{k!} C\textsubscript{k}}}$ \par 
where  C\textsubscript{1} = $<\Delta V>$,  C\textsubscript{2} = $<\Delta V \textsuperscript{2}>$ -  $<\Delta V>$\textsuperscript{2}, and  C\textsubscript{3} = $<\Delta V \textsuperscript{3}>$  - 3 $<\Delta V \textsuperscript{2}>$  + 2 $<\Delta V>$\textsuperscript{3}
\item Free Energy, F(A) =  - k\textsubscript{B} T ln(p(A))
\end{itemize}
```
