%%                           fxn_compdotvar.m
% Alistair Boettiger                                Date Begun: 02/14/11
%                                                   Last Modified: 02/14/11

%% Called by:
% anlz_snail_dots, 02/14/11
% 

function [mRNA1_plot,mRNA2_plot,var1,var2] = fxn_compdotvar(NucLabeled,conn_map,mRNA_sadj1,mRNA_sadj2, Nnucs)

    
    [h,w] = size(NucLabeled);
        mRNA1_plot = zeros(h,w); 
        mRNA2_plot = zeros(h,w);   
        var1 = zeros(1,Nnucs);
        var2 = zeros(1,Nnucs);
        reg_data = regionprops(NucLabeled,'PixelIdxList');
            for k=1:Nnucs
                % Assign all pixels in nucleus N equal to the number of
                % transcripts contained in that nucleus (corrected for
                % area)
                pixes = reg_data(k).PixelIdxList;             
                mRNA1_plot(pixes) = mRNA_sadj1(k);
                mRNA2_plot(pixes) = mRNA_sadj2(k);
                
                % Compute variance among neighbors    k = 30
                Neibs = conn_map(k,:)>1;
                local1_cnts = [mRNA_sadj1(Neibs),mRNA_sadj1(k)];
                var1(k) = std(local1_cnts)/mean(local1_cnts); 
                
                local2_cnts = [mRNA_sadj2(Neibs),mRNA_sadj2(k)];
                var2(k) = std(local2_cnts)/mean(local2_cnts); 
                
%                 C = NucLabeled;  
%                 neib_inds = find(Neibs == 1);
%                 for j = 1:length(neib_inds)
%                     C(C==neib_inds(j)) = 400;
%                 end
%                 figure(1); clf; imagesc(C);
            end
