# ptasensitivity

This repository accompanies the paper:

> A simple optimisation for the MeerKAT Pulsar Timing Array
>
> [Hannah Middleton](https://orcid.org/0000-0001-5532-3622), [Ryan M. Shannon](https://orcid.org/0000-0002-7285-6348), [Matthew Bailes](https://orcid.org/0000-0003-3294-3081), [Andrew D. Cameron](https://orcid.org/0000-0002-2037-4216), [Alessandro Corongiu](https://orcid.org/0000-0002-5924-3141), [Marisa Geyer](https://orcid.org/0000-0002-2822-1919), [Max Jones](https://orcid.org/0009-0005-9124-1348), [Michael Kramer](https://orcid.org/0000-0002-4175-2271), [Matthew T. Miles](https://orcid.org/0000-0002-5455-3474), [Aditya Parthasarathy](https://orcid.org/0000-0002-4140-5616), [Andrea Possenti](https://orcid.org/0000-0001-5902-3731), [Daniel J. Reardon](https://orcid.org/0000-0002-2035-4688).
>
> arXiv: [2505.02524](https://arxiv.org/abs/2505.02524)
> 
> Accepted in MNRAS: [doi.org/10.1093/mnras/staf748](https://doi.org/10.1093/mnras/staf748)

## About this work 

We investigate whether small changes to a pulsar timing array observing schedule can provide gains in signal-to-noise ratio (S/N) for a stochastic gravitational wave background signal from a population of massive black hole binaries. We use the MeerKat Pulsar Timing Array (MPTA) as a test. The approach uses a greedly algorithm to reallocate available integration time between pulsars in the array. The overall observing time dedicated to the MPTA is kept constant so that there is only minimal disruption to the existing observation strategy. We make several assumptions on the gravitational-wave signal and the pulsar noise properties to demonstrate our method. See the paper (linked above) for full details. 

For more information about the MeerKAT Pulsar Timing Array (MPTA) visit [mpta-gw.github.io/](https://mpta-gw.github.io/)

## About this repository

Where to find material: 
 - `data` contains data used for this work, e.g. pulsar positions, pulsar noise assumptions, integration times
 - `snr` scripts for the time-swap analysis 
 - `plotting` scripts used to make plots
 - `paper` material used for the paper including the manuscript and results of the time-swap runs
 - `otherInvestigations` includes a separate investigation used in [Miles et al. (2021)](https://doi.org/10.1093/mnras/stab3549)


## Acknowledgements

The MeerKAT telescope is operated by the South African Radio Astronomy Observatory ([SARAO](https://www.sarao.ac.za/)), which is a facility of the National Research Foundation, an agency of the Department of Science and Innovation. We acknowledge the Wurundjeri People of the Kulin Nation as the Traditional Owners of the land where this work was primarily carried out. Computations were performed on the [OzSTAR](https://supercomputing.swin.edu.au/ozstar) national facility at Swinburne University of Technology. The OzSTAR program receives funding in part from the Astronomy National Collaborative Research Infrastructure Strategy (NCRIS) allocation provided by the Australian Government. HM is grateful to Alberto Vecchio, Paul Brook, Christopher Moore, and Pratyasha Gitika for useful discussions and feedback on the manuscript. This research is supported by the Australian Research Council Centre of Excellence for Gravitational Wave Discovery ([OzGrav](https://www.ozgrav.org/)) (project number  CE170100004), including the OzGrav COVID funding scheme. HM acknowledges support from the UK Space Agency ([UKSA](https://www.gov.uk/government/organisations/uk-space-agency)), Grant No. ST/Y004922/1 and ST/V002813/1 and ST/X002071/1. MK acknowledges significant support from the Max-Planck Society (MPG) and the MPIfR contribution to the PTUSE hardware. AP acknowledges funding from the INAF Large Grant 2022 ``GCjewels'' (P.I. Andrea Possenti) approved with the Presidential Decree 30/2022. This work was supported in part by the ``Italian Ministry of Foreign Affairs and International Cooperation'', grant number ZA23GR03, under the project ``\emph{RADIOMAP- Science and technology pathways to MeerKAT+: the Italian and South African synergy}''. This work has made use of [astropy](https://www.astropy.org/), [numpy](https://numpy.org/), [scipy](https://scipy.org/), [matplotlib](https://matplotlib.org/), [hasasia](https://hasasia.readthedocs.io/en/latest/), and the [Australian Telescope National Facility (ATNF) Pulsar Catalogue](https://www.atnf.csiro.au/research/pulsar/psrcat/).
