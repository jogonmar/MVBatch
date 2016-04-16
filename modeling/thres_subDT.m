function  mp_model2 = thres_subDT(mp_model,Tm)

% Reduction of a MP model in a Down-top fashion.
%
% mp_model2 = thres_subDT(mp_model,Tm) % complete call
%
%
% INPUTS:
%
% mp_model: (structure) input MP model (use the command "help MP_toolbox_h" for
%       more info)
%
% Tm: (1x1) improvement threshold for subdivision.
%
%
% OUTPUTS:
%
% mp_model2: (structure) output MP model (use the command "help MP_toolbox_h" for
%       more info)
%
%
% codified by: José Camacho Páez.
% version: 2.0
% last modification: 09/Feb/2009.

% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;

if ndims(mp_model.arg.xini)~=3, error('Incorrect number of dimensions of xini.'); end;
s = size(mp_model.arg.xini);
if find(s<1), error('Incorrect content of xini.'); end;
if isequal(mp_model.type,'SW-Div'),
    if mp_model.arg.lag<0, error('Incorrect value of lag.'); end;
    if mp_model.arg.lag>s(1)-1, mp_model.arg.lag=s(1)-1; end;
end
if (mp_model.arg.minsize<1||mp_model.arg.minsize>s(1)), error('Incorrect value of minsize.'); end;
if mp_model.arg.n<0, error('Incorrect value of n.'); end;
if (mp_model.arg.prep<0||mp_model.arg.prep>3), error('Incorrect value of prep.'); end;

if ndims(mp_model.phases)~=2, error('Incorrect number of dimensions of phases.'); end;
sp = size(mp_model.phases);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,1:3)<0), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,4:5)<1), error('Incorrect content of phases.'); end;
if find(mp_model.phases(:,4:5)>s(1)), error('Incorrect content of phases.'); end;

if ndims(mp_model.tree)~=2, error('Incorrect number of dimensions of tree.'); end;
sp = size(mp_model.tree);
if (sp(1)<1||sp(2)~=5), error('Incorrect content of tree.'); end;
if find(mp_model.tree(:,1:3)<0), error('Incorrect content of tree.'); end;
if find(mp_model.tree(:,4:5)<1), error('Incorrect content of tree.'); end;
if find(mp_model.tree(:,4:5)>s(1)), error('Incorrect content of tree.'); end;


% Initialization

mp_model2=mp_model;
arg = mp_model.arg;

% Reduction

if Tm <= arg.T, return; end

[cumpress_0,press_0] = crossval3D_s(arg.xini,0,0,mp_model.clu,arg.cross.leave_m,arg.cross.blocks_r,arg.cross.blocks_c,arg.cross.fold_m,arg.prep,arg.cross.order);

mp_model2.tree = thres_sub_recDT(mp_model.tree,Tm,mp_model.arg.absolute,mp_model.arg.gamma,press_0);

[mp_model2.phases, mp_model2.clu] = reduce(mp_model2.tree);

mp_model2.arg.Tm = Tm;

if isequal(mp_model2.type,'SW-Div'),
    mp_model2 = crossvalMP_s(mp_model2);
else
    mp_model2 = crossvalMP_v(mp_model2);
end

end


function tree2 = thres_sub_recDT(tree,Tm,absolute,gamma,val)

cont=true;
tree2=tree;
entries=length(tree(:,4)); 
ind=1;

while cont && ind+1<=length(tree2(:,1)),    
    
    if absolute, 
        baseline = sum(val); 
    else
        baseline = tree2(ind,1); 
    end;
     
    if tree2(ind,5)==tree2(ind+1,5), % Adding a PC
        if ind+1 == length(tree2(:,1)), % Leaf  
            if (tree2(ind,1)-tree2(ind+1,1))/baseline <= Tm,
                tree2 = tree2(1:ind,:);
                ind = ind-1;
                if ind<1, cont=false; end;
            else
                cont=false;
            end            
        else
            ind = ind+1;
        end
        
    else % Adding a phase        
        phases = ind+1;
        end_ind = tree2(ind,5);
        while tree2(phases(end),5) < end_ind,
            new_phase = find(tree2(:,4) == tree2(phases(end),5)+1);
            phases = [phases new_phase(1)];
        end
        
        tree_fp = thres_sub_recDT(tree2(phases(1):phases(2)-1,:),Tm,absolute,gamma,val(1:tree2(phases(2),4)-tree2(1,4)));
        tree3 = [tree_fp];
        for i=2:length(phases)-1,
            tree_fp = thres_sub_recDT(tree2(phases(i):phases(i+1)-1,:),Tm,absolute,gamma,val((tree2(phases(i),4)+1:tree2(phases(i+1),4))-tree2(1,4)));
            tree3 = [tree3;tree_fp];
        end
        tree_sp = thres_sub_recDT(tree2(phases(end):end,:),Tm,absolute,gamma,val(tree2(phases(end),4)-tree2(1,4)+1:end));            
        tree3 = [tree3;tree_sp];
           
        vat2 = sum(tree2(phases,1));
        tree2 = [tree2(1:ind,:);tree3];
        cont = false;
        if length(tree3(:,1))==length(phases), % Leaf
            
            if gamma*(tree2(ind,1)-vat2)/baseline <= Tm,
                tree2 = tree2(1:ind,:);
                ind = ind-1;
                if ind>0, cont=true; end;
            else
                cont=false;
            end
        end

    end
end        

end
