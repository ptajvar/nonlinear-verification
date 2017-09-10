classdef MILPC
    properties
        types
        Variables  %% nV*4 cell matrix - COL1: name COL2: type COL3: length (for arrays) COL4: higher and lower bound
        Constraints %% nC*4 cell matrix - COL1: variables COL2: factors COL3: right hand-side value COL4: Variable number (for arrays)
        ConstraintsEQ
        CostFunction %% 1*3 cell matrix - COL1: variables COL2: factors COL3: Variable number (for arrays)
        nVariables %% number of variables
        nConstraints %% number of constraints
        nConstraintsEQ
        map %% nV*2 cell matrix - COL1: name COL2: its position in MILP
    end
    methods
        function obj = MILPC()
            obj.types = {'double','int','bin'};
            obj.Variables = struct([]);
            obj.Constraints = struct([]);
            obj.ConstraintsEQ = struct([]);
            obj.nVariables = 0;
            obj.nConstraints = 0;
            obj.nConstraintsEQ = 0;
            obj.map = struct([]);
        end
        function obj = Add_Variable(obj, name, type, number, bounds)
            if nargin<5
                bounds = [-inf,inf];
            end
            if nargin < 4
                number = 1;
                if nargin < 3
                    error('define variable type: ''double'', ''int'' or ''bin''');
                end
            end
            knownType = false;
            for i = 1:length(obj.types)
                if strcmp(obj.types{i},type)
                    knownType = true;
                    break;
                end
            end
            if ~knownType
                error('Unknown type: Use ''double'', ''int'' or ''bin''');
            end
            if obj.hasName(name)
                error('Name already exists');
            end
            obj.Variables(end+1).name = name;
            obj.Variables(end).type = type;
            obj.Variables(end).nMembers = number;
            if strcmp(type,'bin')
                obj.Variables(end).bound = repmat([0,1],number,1);
            else
                S = size(bounds);
                if S(1)==1
                    obj.Variables(end).bound = repmat(bounds,number,1);
                elseif S(1)==number
                    obj.Variables(end).bound = bounds;
                else
                    error('bound matrix size should be either (1*2) or (arrayLength*2)');
                end
            end
            obj.map(end+1).name = name;
            obj.map(end).number = obj.nVariables+1;
            obj.nVariables = obj.nVariables+number;
        end
        
        function obj = Add_Constraint(obj, names, factors, rightHand)
            for i = 1:length(names)
                name = names{i};
                if ~iscell(name)
                    name = {name,1}; %#ok<*AGROW>
                    names{i} = name;
                end
                if ~obj.hasName(name{1})
                    err_string = ['Variable name ', name{1}, 'does not exist.'];
                    error(err_string);
                end
            end
            nVar = 0;
            for i =1:length(names)
                if length(names{i})==2
                    nVar = nVar+length(names{i}{2});
                else
                    nVar = nVar+1;
                end
            end
            if nVar~=length(factors)
                error('Variable names and factors length mismatch');
            end
            obj.nConstraints = obj.nConstraints+1;
            obj.Constraints(obj.nConstraints).names = names;
            obj.Constraints(obj.nConstraints).factors = factors;
            obj.Constraints(obj.nConstraints).rightHand = rightHand;
        end
        
        function obj = Add_ConstraintEQ(obj, names, factors, rightHand)
            for i = 1:length(names)
                name = names{i};
                if ~iscell(name)
                    name = {name,1}; %#ok<*AGROW>
                    names{i} = name;
                end
                if ~obj.hasName(name{1})
                    err_string = ['Variable name ', name{1}, 'does not exist.'];
                    error(err_string);
                end
            end
            nVar = 0;
            for i =1:length(names)
                if length(names{i})==2
                    nVar = nVar+length(names{i}{2});
                else
                    nVar = nVar+1;
                end
            end
            if nVar~=length(factors)
                error('Variable names and factors length mismatch');
            end
            obj.nConstraintsEQ = obj.nConstraintsEQ+1;
            obj.ConstraintsEQ(obj.nConstraintsEQ).names = names;
            obj.ConstraintsEQ(obj.nConstraintsEQ).factors = factors;
            obj.ConstraintsEQ(obj.nConstraintsEQ).rightHand = rightHand;
        end
        
        function obj = Set_Cost(obj, names, factors)
            for i = 1:length(names)
                name = names{i};
                if ~iscell(name)
                    name = {name,1};
                    names{i} = name;
                end
                if ~obj.hasName(name{1})
                    err_string = ['Variable name ', name{1}, 'does not exist.'];
                    error(err_string);
                end
            end
            nVar = 0;
            for i =1:length(names)
                if length(names{i})==2
                    nVar = nVar+length(names{i}{2});
                else
                    nVar = nVar+1;
                end
            end
            if nVar~=length(factors)
                error('Variable names and factors length mismatch');
            end
            obj.CostFunction.names = names;
            obj.CostFunction.factors = factors;
        end
        
        function out = getMILP(obj)
            f = zeros(1,obj.nVariables);
            cnt = 1;
            for i = 1:length(obj.CostFunction.names)
                place = obj.mapFind(obj.CostFunction.names{i}{1})+...
                    obj.CostFunction.names{i}{2}-1;
                f(place) = obj.CostFunction.factors(cnt:cnt+length(place)-1);
                cnt = cnt+length(place);
            end
            intcon = zeros(1,obj.nVariables);
            for i = 1:length(obj.Variables)
                if ~strcmp(obj.Variables(i).type,'double')
                    start = obj.mapFind(obj.Variables(i).name);
                    place = start:start+obj.Variables(i).nMembers-1;
                    intcon(place) = 1;
                end
            end
            intcon = find(intcon);
            A = zeros(obj.nConstraints,obj.nVariables);
            b = zeros(obj.nConstraints,1);
            
            for i = 1:obj.nConstraints
                names = obj.Constraints(i).names;
                fact_cnt = 1;
                for j = 1:length(names)
                    base = obj.mapFind(obj.Constraints(i).names{j}{1});
                    for k = obj.Constraints(i).names{j}{2}
                        place =  base + k - 1;
                        A(i,place) = obj.Constraints(i).factors(fact_cnt);
                        fact_cnt = fact_cnt + 1;
                    end
                end
                b(i) = obj.Constraints(i).rightHand;
            end
            
            Aeq = zeros(obj.nConstraintsEQ,obj.nVariables);
            beq = zeros(obj.nConstraintsEQ,1);
            for i = 1:obj.nConstraintsEQ
                names = obj.ConstraintsEQ(i).names;
                fact_cnt = 1;
                for j = 1:length(names)
                    base = obj.mapFind(obj.ConstraintsEQ(i).names{j}{1});
                    for k = obj.ConstraintsEQ(i).names{j}{2}
                        place =  base + k - 1;
                        Aeq(i,place) = obj.ConstraintsEQ(i).factors(fact_cnt);
                        fact_cnt = fact_cnt + 1;
                    end
                end
                beq(i) = obj.ConstraintsEQ(i).rightHand;
            end
            
            lb = zeros(obj.nVariables,1);
            ub = zeros(obj.nVariables,1);
            for i = 1:length(obj.Variables)
                start = obj.mapFind(obj.Variables(i).name);
                finish = obj.mapFind(obj.Variables(i).name)+...
                    obj.Variables(i).nMembers - 1;
                lb(start:finish) = obj.Variables(i).bound(:,1);
                ub(start:finish) = obj.Variables(i).bound(:,2);
            end
            options = optimoptions('intlinprog','Display','final');
            options.TolInteger = 1e-6;
            out = struct('f',f,'intcon',intcon,...
                'Aineq',A,'bineq',b,...
                'Aeq',Aeq,'beq',beq,...
                'lb',lb,'ub',ub,'options',options,...
                'solver','intlinprog');
        end
        
        function solution = getSol(obj)
            prob = obj.getMILP;
            prob.ctype = char(ones(1,length(prob.lb))*'C');
            prob.ctype(prob.intcon) = 'B';
%             sol = intlinprog(prob);
            sol = cplexmilp(prob);
            solution = cell(obj.nVariables,2);
            num = 1;
            for i = 1:length(obj.Variables)
                for j = 1:obj.Variables(i).nMembers
                    solution{num,1} = [obj.Variables(i).name ' ' num2str(j)];
                    solution{num,2} = sol(num);
                    num = num+1;
                end
            end
            
        end
        
        function res = hasName(obj,name)
            res = false;
            for i = 1:length(obj.Variables)
                if strcmp(obj.Variables(i).name,name)
                    res = true;
                end
            end
        end
        function place = mapFind(obj,name)
            for i = 1:length(obj.Variables)
                if strcmp(name,obj.map(i).name)
                    place = obj.map(i).number;
                    return;
                end
            end
            place = -1;
        end
    end
end