# PairedStim ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.04.490684v1))

Analysis for the paired stimulation for the induction of spike-timing dependent plasticty (STDP). Although STDP via electrical stimulation has been demonstrated both *in vivo* and *in vitro*, it often requires a lot of stimulation and is inconsistent for different tested channel pairs. One component of electrical stimulation that is largely overlooked in these studies is the inhibitory response. As a result, this study expanded on previous literature in two major ways: 
1. Deliver paired stimulation with much smaller temporal resolution. The single unit response to stimulation has dynamics that change at sub-millisecond scales, but most STDP experiments use static inter-stimulus intervals. We tested 0.1 ms to 20 ms in 0.5 ms steps. 
2. Use the single unit responses as a measure of connectivity. Most previous literature use macroscopic measures such as stimulus induced movement or most commonly cortico-cortical evoked potentials (CCEPs). However, movements induced by stimulation require the activation of the entire corticospinal tract, and we still don't have a complete understanding of the complex, multiphasic CCEPs. 

## Analyses Performed
- Create a session list (SL) of all experiments that contains sorted spikes, evoked spikes, inhibitory response, and CCEPs. 
- Assess changes in responses (excitatory, inhibitory, CCEPs) with respect to the ISI.
  - I also introduce "normalized ISIs" which are ISIs divided by the timing of the excitatory response.
- Assess changes in control experiments.
- Analyze how each measure changes over time for each experimental condition to uncover underlying mechanisms. 
