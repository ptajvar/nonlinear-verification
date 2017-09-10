function writePWA(pwa)

nu=pwa{1}.nu; % number of inputs
nx=pwa{1}.nx; % number of states
ny=pwa{1}.ny; % number of outputs
% nd=pwa{1}.numHyp;
umin=pwa{1}.umin;
umax=pwa{1}.umax;
xmin=pwa{1}.xmin;
xmax=pwa{1}.xmax;

pwaFileName = [pwa{1}.filename,'_nonlin','.pwa'];
fid = fopen(pwaFileName,'w+');



% ------------------ INPUT -------------------------
fprintf(fid, 'INPUT{');
for ii=1:nu
    fprintf(fid, sprintf('\n   REAL u%d[%.2f,%.2f];',ii,umin(ii),umax(ii)));
end
for ii=1:nx
    fprintf(fid, sprintf('\n   REAL w%d[%s,%s];',ii,num2str(-1),num2str(1)));
end
fprintf(fid, '\nINPUT}');
% ------------------ STATE -------------------------
fprintf(fid, '\nSTATE{');
for ii=1:nx
    fprintf(fid, sprintf('\n   REAL x%d[%.2f,%.2f];',ii,xmin(ii),xmax(ii)));
end
fprintf(fid, '\nSTATE}');
% ------------------ OUTPUT ------------------------
fprintf(fid, '\nOUTPUT{');
for ii=1:nx
    fprintf(fid, sprintf('\n   REAL y%d;',ii));
end
fprintf(fid, '\nOUTPUT}');
% ------------------ DYNAMICS -----------------------
fprintf(fid, '\nDYNAMIC{');
for ii=1:length(pwa)
    if pwa{ii}.repeatedA~=1
        fprintf(fid, '\n');
        fprintf(fid, [pwa{ii}.Aindex, '\n']);
        strA='[';
        for rowx=1:nx
            strA = [strA '['];
            for colx=1:nx
                strA = [strA sprintf(' %.3f',pwa{ii}.A(rowx,colx))];
            end
            strA = [strA ' ]'];
        end
        fprintf(fid, sprintf('%s]',strA));
    end
end
I = eye(nx);
for ii=1:length(pwa)
    if pwa{ii}.repeatedB~=1
        fprintf(fid, '\n');
        fprintf(fid, [pwa{ii}.Bindex '\n']);
        strB='[';
        for rowx=1:nx
            strB = [strB '['];
            for colx=1:nx+nu+1
                strB = [strB sprintf(' %.3f',pwa{ii}.B(rowx,colx))];
            end
            strB = [strB ' ]'];
        end
        fprintf(fid, sprintf('%s]',strB));
    end
end
fprintf(fid, '\nC1\n');
strC = '[';
for ii=1:ny
    strC = [strC '['];
    for jj=1:nx
        strC = [strC sprintf(' %.3f',pwa{ii}.C(ii,jj))];
    end
    strC = [strC ' ]'];
end
fprintf(fid, sprintf('%s]',strC));

fprintf(fid, '\nD1\n');
strD = '[';
for ii=1:ny
    strD = [strD '['];
    for jj= 1:nu+nx+1
        strD = [strD sprintf(' %.3f',pwa{ii}.D(ii,jj))];
    end
    strD = [strD ' ]'];
end
fprintf(fid, sprintf('%s]',strD));
fprintf(fid, '\nDYNAMIC}');

% ---------------------------   RELATIONAL ---------------------------%
fprintf(fid, '\nRELATIONAL{\n');
% for ii=1:nx
%     for jj=1:length(pwa{1}.grid{ii})-2
%         strHyp = sprintf('1*x%d<%.4f\n',ii,pwa{1}.grid{ii}(end-jj));
%         fprintf(fid, strHyp);
%     end
% end
% for ii=1:nu
%     for jj=1:length(pwa{1}.grid{nx+ii})-2
%         strHyp = sprintf('1*u%d<%.4f\n',ii,pwa{1}.grid{nx+ii}(end-jj)) ;
%         fprintf(fid, strHyp);
%     end
% end
fprintf(fid, 'RELATIONAL}');

% ---------------------------   RELATIONAL ---------------------------%
fprintf(fid, '\nREGION{\n');
% for ii=1:length(pwa)
%     strDelta = '[';
%     for jj=1:length(pwa{ii}.delta)
%         strDelta = [strDelta sprintf(' %d', pwa{ii}.delta(jj))];
%     end
%     strDelta = [strDelta ' ]\n'];
%     fprintf(fid, strDelta);
%     fprintf(fid, [pwa{ii}.matrices '\n']);
% end

fprintf(fid, 'REGION}');
fclose all;
%     copyfile(pwaFileName,'C:\Users\Riccardo\Dropbox\Toolbox\C_Code\models');
end