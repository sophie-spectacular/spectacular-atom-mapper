# spectacular-atom-mapper

This atom-mapper depends on RXNMapper (https://github.com/rxn4chemistry/rxnmapper). The key difference is that we are supporting consistent atom-mapping between multiple reactions with the same set of reactants. This atom-mapping also assigns a mapping to every single atom in the reactants (sometime in the near future, I will implement a version which does not).

This atom-mapper's intended usage is for graph neural networks, specifically the Weisfeiler-Lehman Network (WLN) model.
