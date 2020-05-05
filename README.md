# Minimal-Reactome
Identify a minimal ‘reactome’ for E. coli iAF1260. A minimal reactome must be such that the removal of ANY reaction from the network results in a growth rate that is lower than 5% of wildtype growth rate.

A minimal ‘reactome’ consists of all the essential reactions such that the removal of any reaction from the network results in rendered cell growth. This network of essential reactions can be obtained by selectively deleting sets of reactions. These sets of reactions can be obtained by running pFBA that minimizes the sum of fluxes across the entire network and classify genes/reactions under certain
conditions as: essential, pFBA optima, Enzymatically Less Efficient (ELE), Metabolically Less
Efficient (MLE) or pFBA no-flux genes.
ZeroFlux reactions are the ones that do not carry any flux in optimal condition and can be easily
deleted. Blocked reactions cannot contribute to the growth of cell and are blocked, thus can be
deleted, without effecting the biomass production. ELE: these reactions require more flux through
enzymatic steps than alternative pathways that meet the same predicted growth rate. Hence these
reactions can be also be deleted to minimize the norm and maximize biomass flux. MLE on the other
hand are reactions requiring a growth rate reduction, if used, and can be deleted as well. 
But if for a reaction one pathway is MLE while the alternative pathway is ELE then the ELE reaction
should be prioritized for deletion, as we seek to minimize the norm. If such a case arises the value of
resultant biomass flux will be equal to zero. Then, ELE can be prioritized and deleted along with
ZeroFlux reactions and blocked reactions. While relevant MLE reactions will be deleted along with
non Jnz reactions.
The minimal reactome size obtained for E. coli iAF1260 model is 401 reactions, in the code given as
“Jminnz_rxn”. Attempt to delete any other reaction from the reactions containing non-zero flux will
result in zero flux or NAN.
