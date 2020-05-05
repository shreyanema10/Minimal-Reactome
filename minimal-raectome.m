initCobraToolbox;
load ("iAF1260.mat")
model = iAF1260;
% making a copy of model
WT_Model = model;

% defining a growth rate that is lower than 5% of wild-type growthrate. 
Growthrate = 0.95;
WT_FBA = optimizeCbModel(model,'max','one');
WT_f = WT_FBA.f
% changing biomass flux to 5% lower than wild-type flux
model.lb(find(model.c)) = Growthrate*WT_FBA.f;

% running pFBA that minimizes the sum of fluxes across the entire network
% and classify genes and reactions under certain conditions as: essential, pFBA optima, Enzymatically 
% Less Efficient (ELE), Metabolically Less Efficient (MLE) or pFBA no-flux genes


[GeneClasses RxnClasses modelIrrevFM] = pFBA(model, 'geneoption',0, 'tol',1e-7);
% ZeroFlux reactions are the ones that do not carry any flux in optimal condition and can be easily deleted
% Blocked reactions can not contribute to the growth of cell and are
% blocked can be very well deleted
% ELE: these reactions require more flux through enzymatic steps than alternative pathways 
% that meet the same predicted growth rate. Hence these reaction can be
% deleted to minimize the norm and maximize biomass flux.
% MLE: these reactions requiring a growth rate reduction if used and can be
% deleted. 
RxnsToDelete = [RxnClasses.ZeroFlux_Rxns;RxnClasses.Blocked_Rxns; RxnClasses.MLE_Rxns; RxnClasses.ELE_Rxns];
model = changeRxnBounds(model, RxnsToDelete, 0, 'l');
model = changeRxnBounds(model, RxnsToDelete, 0, 'u');
new_FBA = optimizeCbModel(model);
% But if for a reaction one pathway is MLE while the other is ELE
% then the ELE reaction should be prioritized for deletion. If such case
% arises the value of new_FBA.f will be equal to zero.

if new_FBA.f< Growthrate*WT_FBA.f
    % ZeroFlux reaction, Blocked reactions and ELE can be deleted 
    % leaving MLE
    RxnsToDelete = [RxnClasses.ZeroFlux_Rxns;RxnClasses.Blocked_Rxns; RxnClasses.ELE_Rxns];
    model = changeRxnBounds(model, RxnsToDelete, 0, 'l'); 
    model = changeRxnBounds(model, RxnsToDelete, 0, 'u');
    new_FBA = optimizeCbModel(model);
end

count = 0;
% defining flux tolerence value
tol = 1e-07;
mat_Jnz = []; % matrix to store Jnz values
JnzRxnscount = []; % array to store minimal reactome reaction count
for i = 1:length(RxnClasses.pFBAOpt_Rxns)
    newmod = model;
    newmod = changeRxnBounds(model, RxnClasses.pFBAOpt_Rxns(i,1), 0, 'l');
    newmod = changeRxnBounds(model, RxnClasses.pFBAOpt_Rxns(i,1), 0, 'u');
    solution = optimizeCbModel(newmod, 'max', 'zero');
    if (solution.f > 0) && solution.stat == 1
        % finding Jnz reactions for each pFBA optimal reaction 
        Jnz = abs(solution.x) > tol;
        % delete non Jnz reactions
        newmod = changeRxnBounds(model, model.rxns(~Jnz), 0, 'l');
        newmod = changeRxnBounds(model, model.rxns(~Jnz), 0, 'u');
        count = count+1
        newsol = optimizeCbModel(newmod,'max','zero');
        if (newsol.f >= Growthrate*WT_FBA.f) && newsol.stat == 1
            % appending only after feasibility check
            mat_Jnz(size(mat_Jnz,1)+1, :) = (Jnz); 
            JnzRxnscount(length(JnzRxnscount)+1,1) = sum(Jnz);
        end
    end
end
[minreactomesize,I] = min(JnzRxnscount) % minimal reactome 
Jminnz = mat_Jnz(I,:);
% minimal reactome reaction array
Jminnz_rxn = newmod.rxns(Jminnz'==1);
% minimal reactome
model = changeRxnBounds(model, model.rxns(~Jminnz'), 0, 'l');
model = changeRxnBounds(model, model.rxns(~Jminnz'), 0, 'u');