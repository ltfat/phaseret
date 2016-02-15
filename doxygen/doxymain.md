\mainpage

Phase Retrieval Toolbox Library
===============================




General conventions
-------------------

\anchor stft
 \f[
(\mathcal{V}_g f) (\sfreq, \stime)
    = \int_{\mathbb{R}}\! f(\stimearg+\stime)
    \overline{g(\stimearg)} \me^{-\mi 2\pi \sfreq \stimearg  } \,
    \mathrm{d}\stimearg,
    \ \ \sfreq, \stime \in\mathbb{R}\\
\f]

Papers which use this toolbox
\cite ltfatnote044 \cite ltfatnote043 \cite ltfatnote040.

Algorithms
----------

- \ref gla
- \ref rtisila
- \ref legla
- \ref spsi
- \ref heapint

Block processing
----------------


~~~~~~~~~~~~~~~{.c}
int func(int a,int b)
{ 
    return a*b;
}
~~~~~~~~~~~~~~~


> Blockquote
> 2nd line

- More text for this item.
- Further
  + nested list item.
  + another nested item.
- Prd

1. First item.
2. Second item.

This a normal paragraph

    This is a code block

We continue with a normal paragraph again.


Code inline `int prd`

Cross link to \ref stft "(formula)"



