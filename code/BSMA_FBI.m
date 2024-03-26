
function [GlobalMin,X,Fit_store,Time] = BSMA_FBI1(NP,GEN,D,A,trn,vald,TFid)
tic;
%% Initialized population
ObjVal = zeros(1,NP); 
LB=0;
UB=1;
FES=1;
Pop = initialization(NP,D,1,0)>0.5;
for i = 1:NP 
    ObjVal(1,i) = AccSz2(Pop(i,:), A,trn,vald);
end
%% Memorize the best solution
iBest = find(ObjVal==min(ObjVal));
iBest = iBest(end); 
GlobalMin = ObjVal(iBest); 
Xbest = Pop(iBest,:); 
g = 1;
Fit_store = [];
%% Optimization Cycle
PopA = Pop; ObjValA = ObjVal;
PopB = Pop; ObjValB = ObjVal;
PopC = Pop;

Destination_fitness=inf;
AllFitness = inf*ones(NP,1);
weight = ones(NP,D);
z = 0.03;

% f_a = inf*ones(NP,1);
%  f_a(i) = AccSz2(PopA(i,:), A,trn,vald);
while FES <= GEN
    %% Investigation team - team A
    %% Step A1
    f_a = inf*ones(NP,1);
    for i=1:NP %ÐÐ
        Change=fix(rand*D)+1; 
        nb1=floor(rand*NP)+1;      
            while(nb1==i)
                nb1=floor(rand*NP)+1;
            end
        nb2=floor(rand*NP)+1;      
            while(nb1==nb2 || nb2==i)
                nb2=floor(rand*NP)+1;
            end   
        solA=PopA(i,:);
       change=PopA(i,Change)+(PopA(i,Change)-(PopA(nb1,Change)+PopA(nb2,Change))/2)*(rand-0.5)*2; %Eq.(2) in FBI Inspired Meta-Optimization
        solA(Change)= transferFun(solA(Change),change,TFid);
       %        for Change = 1:D
%             if(solA(Change)<=LB(Change))||(solA(Change)>=UB(Change))
%                 solA(Change) = LB(Change) + (UB(Change)-LB(Change))*rand();
%             end
%        end

        if(solA(Change)<=0)||(solA(Change)>=1)
            solA(Change) = 0 + (1-0)*rand();
        end
        

        f_a(i) = AccSz2(solA, A,trn,vald);
%         f_a = AccSz2(solA, A,trn,vald,classifierFhd);
        if f_a(i) <= ObjValA(i)
            PopA(i,:) = solA; 
            ObjValA(i) = f_a(i) ;
            if f_a <= GlobalMin
                Xbest = solA ; 
                GlobalMin = f_a(i); 
            end  
        end
    end
    %% Step A2
    if min(ObjValA) < max(ObjValA)
    prob = probability(ObjValA);
    for i = 1:NP
        if (rand>prob(i))
            r(1) = floor(rand()* NP) + 1;
            while r(1)==i
                r(1) = floor(rand()* NP) + 1;
            end
            r(2) = floor(rand()* NP) + 1;
            while (r(2)==r(1))||(r(2)==i)
                r(2) = floor(rand()* NP) + 1;
            end
             r(3) = floor(rand()* NP) + 1; 
             while (r(3)==r(2))||(r(3)==r(1))||(r(3)==i)
                 r(3) = floor(rand()* NP) + 1; 
             end 
            solA = PopA(i,1:D);
            Rnd = floor(rand()*D) + 1;
            for j = 1:D 
                if (rand()< rand()) || ( Rnd==j)
                    change = Xbest(j) + PopA(r(1),j) + rand() * (PopA(r(2),j) - PopA(r(3),j)); %Eq.(5) in FBI Inspired Meta-Optimization
                    solA(j)= transferFun(solA(j),change,TFid);
                else
                    solA(j) = PopA(i,j);
                end
                
                if(solA(j)<=0)||(solA(j)>=1)
                    solA(j) = 0 + (1-0)*rand();
                end
            end
%             for j = 1:D
%                 if (solA(j)<=LB(j))||(solA(j)>=UB(j))
%                     solA(j) = LB(j) + (UB(j)-LB(j))*rand();
%                 end
%             end

            f_a = AccSz2(solA, A,trn,vald);
            if f_a <= ObjValA(i)
                PopA(i,:) = solA; 
                ObjValA(i) = f_a;
                if f_a <= GlobalMin
                    Xbest = solA ; 
                end  
            end
        end 
    end
    end
    %% Persuing team - team B
    %% Step B1 
    for i = 1:NP 
        SolB = PopB(i,1:D);  
        for j = 1:D
            change = rand()*PopB(i,j) + rand()*(Xbest(j) - PopB(i,j)); %Eq.(6) in FBI Inspired Meta-Optimization
            SolB(j) = transferFun(SolB(j),change,TFid);
            %             if (SolB(j)<LB(j))||(SolB(j)>UB(j))
%                 SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
%             end
            if(SolB(j)<=0)||(SolB(j)>=1)
                SolB(j) = 0 + (1-0)*rand();
            end
        end
        f_b = AccSz2(SolB, A,trn,vald);
        if f_b <= ObjValB(i)
            PopB(i,1:D) = SolB; 
            ObjValB(i) = f_b;
            if f_b <= GlobalMin
                Xbest = SolB ; 
                GlobalMin = f_b; 
            end
        end
    end 
    %% Step B2
    for i = 1:NP 
        rr = floor(rand()* NP) + 1;
        while rr ==i
            rr = floor(rand()* NP) + 1;
        end
        if ObjValB(i) > ObjValB(rr)
            temp = PopB(rr,:) + rand(1,D).* (PopB(rr,:) - PopB(i,:))+ rand(1,D).* (Xbest - PopB(rr,:)); %Eq.(7) in FBI Inspired Meta-Optimization
            for j = 1:D
%                 if (SolB(j)<LB(j))||(SolB(j)>UB(j))
%                     SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
%                 end
                SolB(j) = transferFun(SolB(j),temp(j),TFid);
                if(SolB(j)<=0)||(SolB(j)>=1)
                    SolB(j) = 0 + (1-0)*rand();
                end
            end

        else
            temp = PopB(i,:) + rand(1,D).* (PopB(i,:)-PopB(rr,:))+ rand(1,D).* (Xbest - PopB(i,:)); %Eq.(8) in FBI Inspired Meta-Optimization
            for j = 1:D
%                 if (SolB(j)<LB(j))||(SolB(j)>UB(j))
%                     SolB(j) = LB(j) + (UB(j)-LB(j))*rand();
%                 end
                SolB(j) = transferFun(SolB(j),temp(j),TFid);
                if(SolB(j)<=0)||(SolB(j)>=1)
                    SolB(j) = 0 + (1-0)*rand();
                end
            end

        end
        f_b = AccSz2(SolB, A,trn,vald);
        if f_b <= ObjValB(i)
            PopB(i,:) = SolB; 
            ObjValB(i) = f_b;
            if f_b <= GlobalMin
                Xbest = SolB ; 
                GlobalMin = f_b; 
            end 
        end
    end 
    
     %% Persuing team - team C
     

    for i=1:NP
        AllFitness(i) = AccSz2(PopC(i,:), A,trn,vald);
    end
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(NP);
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;  


    for i=1:NP
        for j=1:D
            if i<=(NP/2)    %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    

    if bestFitness < Destination_fitness
        bestPositions=PopC(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(FES/GEN)+1); %Eq.(2.4)
    b = 1-FES/GEN; 
    
    XX = zeros(1,D);

    for i=1:NP
        if rand<z     %Eq.(2.7)
            XX(1,:) = (UB-LB)*rand+LB;
			for j=1:D
				ss = XX(1,j);
				PopC(i,j)=transferFun(PopC(1,j),ss,TFid);				
			end
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness)); %Eq.(2.2)
            vb = unifrnd(-a,a,1,D);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,D);
            for j=1:D
                r = rand();
                AA = randi([1,NP]);   
                B = randi([1,NP]);
                if r<p          %Eq.(2.1)
                    ss = bestPositions(j)+ vb(j)*(weight(i,j)*PopC(AA,j)-PopC(B,j));
                else
                    ss = vc(j)*PopC(i,j);
                end
                PopC(i,j)=transferFun(PopC(i,j),ss,TFid);
            end
        end
    end
    
    if bestFitness <= GlobalMin
        GlobalMin = bestFitness;
    end
    
%     if rem(g, ishow) == 0
%         fprintf('Generation: %d. Best f: %f.\n', g, GlobalMin);
%     end
    Fit_store(FES) = GlobalMin;
    FES = FES + 1;
end 
%% Result
f = GlobalMin;
X = Xbest;   
Time = toc;
end 

% %% functions
% function f = rescale_matrix(X, LB, UB)
% [NP,D] = size(X);
% f = zeros(NP,D);
% for i = 1:D
%     f(:, i) = LB(i)*ones(NP,1) + (UB(i) - LB(i))*X(:,i);
% end
% end
%% functions
function prob = probability(fObjV)
prob = (max(fObjV)-fObjV)/((max(fObjV)-min(fObjV))); %Eq.(3) in FBI Inspired Meta-Optimization
end







