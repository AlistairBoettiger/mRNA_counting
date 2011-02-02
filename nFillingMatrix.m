

%%
% Alistair Boettiger                                   Date Begun: 01/22/11
% Functionally complete                             Last Modified: 01/22/11


% Build the filling matrix of order n with transition rates that depend on
% both state and direction.  

function [Gu Gl] = nFillingMatrix(n,ks,km)

% n = 5;
% ks  = 1:n+1; % rate constants dependent on n .  % is this right? 
% km = (1:n+1)/sqrt(3); 


Nstates = 2^n;
All_states = dec2bin(0:Nstates-1);  % Nstates x n
G = zeros(Nstates);

for i=1:Nstates  % i = 11
  cur_state = All_states(i,:);
  nbound = length( strfind(cur_state,'1') );
  
  neib_state = cell(1,n); 
  for j=1:n
      % Find the binary and numeric index of all neibhbor states.
     new_bit = num2str( rem(str2double(cur_state(j))+1,2) ); % flip the jth bit  
     neib_state{j} = [cur_state(1:j-1),new_bit,cur_state(j+1:end)]; % concat with the other bits
     
      dir = sum( neib_state{j} - cur_state) ; % this is +1 for up and -1 for a down transition
     % convert binary string into the index for that neighboring state. 
     ind_neib = bin2dec(neib_state{j})+1;  
     
     if dir == 1
     G(i,ind_neib) = ks(nbound+1);  
     elseif dir == -1
     G(i,ind_neib) = km(nbound+1);  
     else
         disp('error');
     end
  end

end

Gu = triu(G);  Gl = tril(G); 

gsu = sum(Gu,2);  Gu = Gu-diag(gsu); 
gsl = sum(Gl,2);  Gl = Gl-diag(gsl);

figure(2); clf; imagesc(Gu+Gl);



%   %% Compute the connectivity for the filling matrix of order n. 
% 
% n = 5;
% Nstates = 2^n;
% All_states = dec2bin(0:Nstates-1);  % Nstates x n
% 
% 
% G = zeros(Nstates);
% 
% for i=1:Nstates
%   cur_state = All_states(i,:);
%   
%   neib_state = cell(1,n); 
%   for j=1:n
%       % Find the binary and numeric index of all neibhbor states.
%       new_bit = num2str( rem(str2double(cur_state(j))+1,2) ); % flip the jth bit  
%      neib_state{j} = [cur_state(1:j-1),new_bit,cur_state(j+1:end)]; % concat with the other bits
%      
%      % convert binary string into the index for that neighboring state. 
%      ind_neib = bin2dec(neib_state{j})+1;  
%      G(i,ind_neib) = 1;  
%   end
% 
% end
% 
% figure(1); clf; imagesc(G);