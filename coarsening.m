function [Sol,eta] = coarsening(Sol,eta,theta)
% COARSENING coarsens the current mesh
%
% USAGE
%    [Sol,eta] = coarsening(Sol,eta,theta)
%
% INPUT 
%      Sol:  current mesh data
%      eta:  error indicator for each triangle
%    theta:  parameter in (0,1). 
%            We mark minimal number of triangles M such that
%               \sum_{T \in M} \eta_T < \theta*\sum\eta_T
%
% OUTPUT
%      Sol:  new mesh data after coarsening
%      eta:  new error indicator for each triangle
%
% NOTE
%    Only one level coarsening
%
% REFERENCE
%    AFEM@matlab Manual --> Algorithms --> Coarsening     
%

% L. Chen & C. Zhang 11-15-2006
% Modified by Vaibhav Joshi, Rajeev K. Jaiman

%--------------------------------------------------------------------------
% Construct data structure
%--------------------------------------------------------------------------
N = size(Sol.node,1); 
NT = size(Sol.elem,1);
dualEdge = sparse(Sol.elem(:,[1,2,3]),Sol.elem(:,[2,3,1]),[1:NT,1:NT,1:NT]);
valence = accumarray(Sol.elem(:),ones(3*NT,1),[N 1]);
eta_p = accumarray(Sol.elem(:),[eta;eta;eta],[N 1]);
 
%--------------------------------------------------------------------------
% Mark nodes for coarsen
%--------------------------------------------------------------------------
total = sum(eta_p); maxeta = max(eta_p);
current = 0; marker = uint8(zeros(N,1)); 
newestNode = unique(Sol.elem(:,1)); % newest vertices of all elements
coarseNode = find(Sol.type==1); % nodes in the initial triangulation
newestNode = setdiff(newestNode,coarseNode);
isGood = (valence(newestNode) == 2) | (valence(newestNode) == 4);
goodNode = newestNode(isGood);
% A node is called 'good' if it can be removed without breaking the mesh
% comformity. A node is 'good' if and only if
%  0. it is not a node of marco elements in the initial triangluation;
%  1. it is the newest node in each of the elements with it as a node;
%  2. if it is an interior node, the valence is 4;
%  3. if it is on the boundary, the valence is 2.
[temp,ix] = sort(eta_p(goodNode));
for i = 1:length(goodNode)
    p = goodNode(ix(i));
    if (eta_p(p) > 0.2*maxeta), break; end % only coarsen where error small
    if (current + eta_p(p) > theta*total), break; end
    marker(p) = 1;    % mark p to be removed
    Sol.type(p) = 0; % node p is non active now
    current = current + eta_p(p);
end
if (sum(marker)==0), return; end % no node can be coarsened

%--------------------------------------------------------------------------
% Coarsen the mesh by removing marked nodes
%--------------------------------------------------------------------------
for t = 1:NT
    if (Sol.elem(t,1)>0) % otherwise t is already deleted
        p = Sol.elem(t,1);
        if (marker(p)==1) % p is a marked good point
            brother = dualEdge(Sol.elem(t,2),p);
            Sol.elem(t,1) = Sol.elem(t,2);
            Sol.elem(t,2) = Sol.elem(t,3);
            Sol.elem(t,3) = Sol.elem(brother,2);
            Sol.elem(brother,1) = 0; % brother is going to be deleted
            eta(t) = eta(t) + eta(brother); % new eta = eta_t + eta_brother
        end % endif for all markd good nodes
    end
end % end of for loop on all elements

%--------------------------------------------------------------------------
% Clear element and eta arrarys
%--------------------------------------------------------------------------
ix = (Sol.elem(:,1) == 0); 
Sol.elem(ix,:) = []; 
eta(ix,:) = [];

%--------------------------------------------------------------------------
% Update boundary edges
%--------------------------------------------------------------------------
markedNode = goodNode((marker(goodNode)==1)&(valence(goodNode)==2));
Sol.Dirichlet = updatebd(Sol.Dirichlet, markedNode);
Sol.Neumann   = updatebd(Sol.Neumann, markedNode);

%--------------------------------------------------------------------------
% End of function COARSENING
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sub functions called by COARSENING
%--------------------------------------------------------------------------
function bdEdge = updatebd(bdEdge, markedNode)
% UPDATEDBD coarsen the boundary edges
%
% USAGE
%    bdEdge = updatebd(bdEdge, markedNode)
%
% INPUT
%     bdEdge:  old boundary edges
% markedNode:  node marked for coarsening
%
% OUTPUT
%     bdEdge:  new boundary edges after coarsening
%

if (isempty(bdEdge) || isempty(markedNode)), return; end
NB = size(bdEdge,1); 
neig = sparse(bdEdge,ones(NB,1)*[2,1],[1:NB,1:NB]);
% skeleton of edges around a node
for k = 1: length(markedNode)
    i = neig(markedNode(k),1);
    j = neig(markedNode(k),2);
    % o --- i --- o --- j --- o
    %             ^
    %         noActive(k)
    bdEdge(i,2) = bdEdge(j,2);
    bdEdge(j,1) = 0;
end
ix = (bdEdge(:,1)==0); 
bdEdge(ix,:) = [];
%--------------------------------------------------------------------------
% End of function UPDATEDBD
%--------------------------------------------------------------------------
