function [assignment,cost] = munkres(costMat)
%匈牙利算法解决线性分配问题
%
%[ASSIGN,COST] = munkres(COSTMAT)返回一个最优的列索引。
%ASSIGN表示匹配行索引，COST表示基于分配的最小花费值。

assignment = zeros(1,size(costMat,1));
cost = 0;

validMat = costMat == costMat & costMat < Inf;
bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
costMat(~validMat) = bigM;

% costMat(costMat~=costMat)=Inf;
% validMat = costMat<Inf;
validCol = any(validMat,1);
validRow = any(validMat,2);

nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return
end

maxv=10*max(costMat(validMat));

dMat = zeros(n) + maxv;
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);

%*************************************************
% 算法开始
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   第一步:从每一行中减去最小行数。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minR = min(dMat,[],2);
minC = min(bsxfun(@minus, dMat, minR));

%**************************************************************************  
%   步骤2：在dMat中查找0。如果其列或行中没有标记的零，则从零开始。对每个0重复
%**************************************************************************
zP = dMat == bsxfun(@plus, minC, minR);

starZ = zeros(n,1);
while any(zP(:))
    [r,c]=find(zP,1);
    starZ(r)=c;
    zP(r,:)=false;
    zP(:,c)=false;
end

while 1
%**************************************************************************
%  步骤3:用星号为零的符号覆盖每一列。如果所有列都是然后覆盖匹配是最大的
%**************************************************************************
    if all(starZ>0)
        break
    end
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;
    coverRow = false(n,1);
    primeZ = zeros(n,1);
    [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
    while 1
        %**************************************************************************
        %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
        %           zero in the row containing this primed zero, Go to Step 5.  
        %           Otherwise, cover this row and uncover the column containing 
        %           the starred zero. Continue in this manner until there are no 
        %           uncovered zeros left. Save the smallest uncovered value and 
        %           Go to Step 6.
        %**************************************************************************
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);
            uZc = cIdx(1);
            primeZ(uZr) = uZc;
            stz = starZ(uZr);
            if ~stz
                Step = 5;
                break;
            end
            coverRow(uZr) = true;
            coverColumn(stz) = false;
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            % *************************************************************************
            % 步骤6:将最小未覆盖值添加到每一行的每个元素中，并从每一列的每个元素中减去它。回到步骤4，不要改变任何星星、质数或覆盖线。
            %**************************************************************************
            [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
            minC(~coverColumn) = minC(~coverColumn) + minval;
            minR(coverRow) = minR(coverRow) - minval;
        else
            break
        end
    end
    %**************************************************************************
    % 步骤5:构造一系列交替的素数零和星号零，如下所示:让Z0表示步骤4中发现的未发现的素数零。设Z1表示Z0列中的星号0(如果有的话)。让Z2表示Z1行中的素数0(总是有1)。把这个级数的每一颗星都去星，把这个级数的每一撇都去星，擦掉所有的质数，然后找出矩阵中的每一行。回到步骤3。
    %**************************************************************************
    rowZ1 = find(starZ==uZc);
    starZ(uZr)=uZc;
    while rowZ1>0
        starZ(rowZ1)=0;
        uZc = primeZ(rowZ1);
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);
        starZ(uZr)=uZc;
    end
end

% 计算分配的总开销
rowIdx = find(validRow);
colIdx = find(validCol);
starZ = starZ(1:nRows);
vIdx = starZ <= nCols;
assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
pass = assignment(assignment>0);
pass(~diag(validMat(assignment>0,pass))) = 0;
assignment(assignment>0) = pass;
cost = trace(costMat(assignment>0,assignment(assignment>0)));

function [minval,rIdx,cIdx]=outerplus(M,x,y)
ny=size(M,2);
minval=inf;
for c=1:ny
    M(:,c)=M(:,c)-(x+y(c));
    minval = min(minval,min(M(:,c)));
end
[rIdx,cIdx]=find(M==minval);
