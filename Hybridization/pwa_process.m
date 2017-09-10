function pwa_post = pwa_process(pwa_pre,filename) % post-processing pwa to be ready for 'writePWA'
if nargin<2
    filename = char(datetime);
    filename(filename=='-' | filename==':') = '_';
end
pwa_post = pwa_pre;
pwa_post{1}.filename = filename;
cA = 0;
cB = 0;
for ii=1:length(pwa_post)
    pwa_post{ii}.Aindex=sprintf('A%d',ii-cA);
    for jj=1:ii-1
        if isequal(pwa_post{ii}.A,pwa_post{jj}.A) % fill the matricesA with all the available A matrices (even those repeated)
            pwa_post{ii}.Aindex = pwa_post{jj}.Aindex;
            pwa_post{ii}.repeatedA = 1;
            cA=cA+1;
            break;
        end
    end
    pwa_post{ii}.Bindex=sprintf('B%d',ii-cB);
    for jj=1:ii-1
        if isequal(pwa_post{ii}.B,pwa_post{jj}.B)
            pwa_post{ii}.Bindex = pwa_post{jj}.Bindex;
            pwa_post{ii}.repeatedB = 1;
            cB=cB+1;
            break;
        end
    end
    pwa_post{ii}.matrices=[pwa_post{ii}.Aindex,',',pwa_post{ii}.Bindex,',C1,D1'];
end
% numA = prod(nbin)-cA;
% numB = prod(nbin)-cB;