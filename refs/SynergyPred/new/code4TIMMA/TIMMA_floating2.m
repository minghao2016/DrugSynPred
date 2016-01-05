% initial list: 1xk set of 0s/1s, where a 1 means the kinase is in the initial set
% protected list: 1xk 0s/1s,1 means kinase will not be removed when added.
% using floating search to select kinases
% NB! as long as the kinase sets are made the same partition of the drugs
% then the LOO err is identical.
% For example: kinase sets (7,8,44,102,211,226,243,251,265,272) make the
% drugs partition the same as kinase sets
% (6,7,8,44,102,211,226,243,251,265,272)

% Tao Xu

function [M,err,k_selected] = TIMMA_floating2(profile,sens, sp, max_k,LOO)
% profile is the drug target profile data
% sens is the drug sensitivity data IC50 ACTarea
% sp: start point of the algorithms
% max_k: the upper limit of the number of kinases.

% [drugs,k] = size(S);
tic();
% use leave one out prediction err as the selection criterion
% LOO = 1; 
% y = spec';


[drugs,k] = size(profile);
S = profile;
% LOO = 0; 
% k_set = initial_list;
y = sens';
new_initial_list = zeros(1,k);
new_initial_list(sp)=1;

% the criterion function J(X_k)
J = zeros(1,k);
% the set X_k
X = zeros(k,k);

run = 0;
step = 0;
change = 1;
% new_initial_list = initial_list;
global_best_err = 100;

while(change)
    run = run +1;
%   fprintf(1,'Run %d\n',run);
    k_set = new_initial_list;
%   fprintf(1,'Start kinase=%d\n',find(k_set));
    
    while(sum(k_set) <= max_k)
        err = ones(1,k)*inf; % individual significance of the feature
        step = sum(k_set);
        
        k_temp = k_set;
        k_temp2 = (k_temp == 1);
        S_temp = S(:,k_temp2);

        [dummy,err0] = TIMMA2(drugs, S_temp,y,sum(k_temp),LOO);
        
        % LOO of nonzero drugs
		%err0 = err0(find(sum(S_temp,2)~=0));
        
        J(step) = sum(err0); % the score for X_k, lower score is better
        if(isequal(X(:,step),k_set')) % if k_set has been visited
            tmp = find(J==min(J(J>0)));
            k_set = X(:,tmp(1));
            change = 0;
            break
        else
            X(:,step) = k_set; % current X_k
        end
		fprintf(1,'k=%d err=%4.2f\n',step,J(step));
        if (step>1 && J(step)==J(step-1)) % if no improvement
            new_initial_list = zeros(1,k);
            new_initial_list(ind1) = 1;
            if (J(step)<global_best_err)
                global_best_err = J(step);
                change = 1;
				fprintf(1,'Best err: %4.2f\n',global_best_err);
            else 
                change = 0;
            end
            break;
        end
        
        % Step 1: Inclusion
		%[IM_ds,IM_is,IM_os] = search_space(drugs,k_set,S,y);
        [IM_ds,IM_is,IM_os] = search_space(drugs,k_set,S,y);
        for i = 1:k
            % fprintf(1,'\b\b\b%d',i)
            % fprintf_r('%d',i);
            if(~k_set(i))
                % Including new kinase i
                err(i) = 0;
                [err1] = TIMMA_search2(S(:,i),IM_ds,IM_is,IM_os,y,LOO);
				%err1 = err1(find(sum(S(:,[i find(k_set)]),2)~=0)); %LOO error of nonzero drugs
				% the cumulative prediction error of including kinase i
                err(i) = sum(err1); 
            end
        end
        
        % select x_{k+1} the most significant feature with respect to X_k
        [dummy, ind1] = min(err);
        k_set(ind1) = 1; 
		%tmp = k_names0(find(ix2==ind1));
		fprintf(1,'Inclusion kinase=%d err=%4.2f\n',ind1,dummy);
		%fprintf(1,'''%s'',',tmp{:});
		%fprintf(1,'\n');
        
		% Step 2: Conditional exclusion
        % find the least significant feature in the set X_{k+1}
        err = ones(1,k)*inf; % individual significance of removing feature i in X_{k+1}
        for i = 1:k
            k_temp = k_set;
            if(k_set(i)) % if removing i
                err(i) = 0;
                k_temp(i) = 0; % remove kinase i from X_{k+1}
                k_temp2 = (k_temp == 1);
                %S_temp is made of the columns of included
                S_temp = S(:,k_temp2);
                [dummy,err1] = TIMMA2(drugs, S_temp,y,sum(k_temp),LOO);
				%err1 = err1(find(sum(S_temp,2)~=0)); % LOO error of nonzero drugs
                err(i) = sum(err1);
            end
        end
        [worst_err,ind_worst] = min(err);
        if (err(ind1)~=worst_err && worst_err < J(step)) % if x_r is the least significant feature
            k_set(ind_worst) = 0; % exclude x_r
			%tmp = k_names0(find(ix2==ind_worst));
			fprintf(1,'Exclusion kinase=%d err=%4.2f\n',ind_worst,worst_err);
			%fprintf(1,'%s ',tmp{:});
			%fprintf(1,'\n');
            while (sum(k_set)>2) % Step 3: Continuation of conditional exclusion
                err = ones(1,k)*inf;
                for i = 1:k
                    k_temp = k_set;
                    if(k_set(i))
                        err(i) = 0;
                        k_temp(i) = 0;
                        k_temp2 = (k_temp==1);
                        S_temp = S(:,k_temp2);
                        [dummy,err1] = TIMMA2(drugs, S_temp,y,sum(k_temp),LOO);
						% err1 = err1(find(sum(S_temp,2)~=0));% LOO error of nonzero drugs
                        err(i) = sum(err1);
                    end
                end
                [dummy, ind_worst] = min(err);
                if (dummy<J(length(find(k_set))-1)) % if better result found
                    k_set(ind_worst) = 0; 
					% tmp = k_names0(find(ix2==ind_worst));	
					fprintf(1,'Continuing exclusion kinase=%d err=%4.2f\n',ind_worst,dummy);	
					% fprintf(1,'%s ',tmp{:});
                    % fprintf(1,'\n');
                else
                    break
                end                
            end
        end
        % set k = k+1 and return to step 1
    change = 0;    
    end  
end


% the best set
% k_set(ind1) = 0;
temp=find(J==min(J(J>0)));
if(length(temp)>1)

   k_set=X(:,temp(2));
else
   k_set=X(:,temp(1));
end
k_temp = k_set;
k_selected= find(k_set);
k_temp2 = (k_temp == 1);
%S_temp is made of the columns of included
S_temp = S(:,k_temp2);
[M,err] = TIMMA2(drugs, S_temp,y,sum(k_temp),LOO);
toc();


