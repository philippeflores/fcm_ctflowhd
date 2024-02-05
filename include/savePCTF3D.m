function savePCTF3D(varargin)
% This function permits to save the results of PCTF3D for future
% visualizations.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

listBan = {'strVariables','strMethod','dataMarg','dataMargFull','X'};

if size(varargin,2)==0

    listVar = evalin("base","who");

    indKept = ones(size(listVar,1),1);
    for i = 1:size(indKept,1)
        for j = 1:size(listBan,2)
            if strcmpi(listBan{j},listVar{i})==1
                indKept(i) = 0;
            end
        end
    end

    listSave = listVar(indKept==1);

    boolStrFile = evalin("base","exist('strFile','var')");
    if boolStrFile==0
        error("\nDefault filename was chosen but it requires the parameter 'strFile' which does not exist in the current workspace.\n")
    else
        strFile = evalin("base","strFile");
    end

    boolI = evalin("base","exist('I','var')");
    if boolI==0
        error("\nDefault filename was chosen but it requires the parameter 'I' which does not exist in the current workspace.\n")
    else
        I = evalin("base","I");
    end

    boolR = evalin("base","exist('R','var')");
    if boolR==0
        error("\nDefault filename was chosen but it requires the parameter 'R' which does not exist in the current workspace.\n")
    else
        R = evalin("base","R");
    end

    boolCalT = evalin("base","exist('calT','var')");
    if boolCalT==0
        error("\nDefault filename was chosen but it requires the parameter 'calT' which does not exist in the current workspace.\n")
    else
        calT = evalin("base","calT");
    end

    filename = sprintf("./save/%s_I%dR%dT%d.mat",strFile,I,R,size(calT,2));

    if isempty(listSave)==0
        evalin("base",sprintf("save('%s','%s')",filename,listSave{1}))
        for i = 2:size(listSave,1)
            evalin("base",sprintf("save('%s','%s','-append')",filename,listSave{i}))
        end
    else
        error("\nThere are no variables to save.\n")
    end

elseif size(varargin,2)==1 && any(char(varargin{1})=='.')

    arg = char(varargin{1});

    if size(arg,2)>4 && strcmpi(arg(end-3:end),'.mat')==1

        filename = arg;
        flag = 0;
        if size(filename,2)>9 && strcmpi(filename(1:5),'save/')==1
            flag = 1;
        end
        if size(filename,2)>11 && strcmpi(filename(1:7),'./save/')==1
            flag = 1;
        end
        if flag==0
            filename = sprintf("./save/%s",filename);
        end

        listVar = evalin("base","who");

        indKept = ones(size(listVar,1),1);
        for i = 1:size(indKept,1)
            for j = 1:size(listBan,2)
                if strcmpi(listBan{j},listVar{i})==1
                    indKept(i) = 0;
                end
            end
        end

        listSave = listVar(indKept==1);

        if isempty(listSave)==0
            evalin("base",sprintf("save('%s','%s')",filename,listSave{1}))
            for i = 2:size(listSave,1)
                evalin("base",sprintf("save('%s','%s','-append')",filename,listSave{i}))
            end
        else
            error("\nThere are no variables to save.\n")
        end

    else
        error("\nThe input argument '%s' is not valid. It should end with '.mat'.\n",arg)
    end
else

    flagFile = 0;
    for i = 1:size(varargin,2)
        varargin{i} = char(varargin{i});
        if flagFile==0 && size(varargin{i},2)>4 && strcmpi(varargin{i}(end-3:end),'.mat')
            flagFile = i;
        end
    end

    if flagFile>0
        filename = varargin{flagFile};
        flag = 0;
        if size(filename,2)>9 && strcmpi(filename(1:5),'save/')==1
            flag = 1;
        end
        if size(filename,2)>11 && strcmpi(filename(1:7),'./save/')==1
            flag = 1;
        end
        if flag==0
            filename = sprintf("./save/%s",filename);
        end
    else
        boolStrFile = evalin("base","exist('strFile','var')");
        if boolStrFile==0
            error("\nDefault filename was chosen but it requires the parameter 'strFile' which does not exist in the current workspace.\n")
        else
            strFile = evalin("base","strFile");
        end

        boolI = evalin("base","exist('I','var')");
        if boolI==0
            error("\nDefault filename was chosen but it requires the parameter 'I' which does not exist in the current workspace.\n")
        else
            I = evalin("base","I");
        end

        boolR = evalin("base","exist('R','var')");
        if boolR==0
            error("\nDefault filename was chosen but it requires the parameter 'R' which does not exist in the current workspace.\n")
        else
            R = evalin("base","R");
        end

        boolCalT = evalin("base","exist('calT','var')");
        if boolCalT==0
            error("\nDefault filename was chosen but it requires the parameter 'calT' which does not exist in the current workspace.\n")
        else
            calT = evalin("base","calT");
        end

        filename = sprintf("./save/%s_I%dR%dT%d.mat",strFile,I,R,size(calT,2));
    end

    for i = 1:size(varargin,2)
        if i~=flagFile
            if strcmpi(varargin{i},'data')
                indKept = ones(size(listBan,2),1);
                for j = 1:size(listBan,2)
                    if strcmpi(listBan{j},'X')
                        indKept(j) = 0;
                    end
                end
                listBan = listBan(indKept==1);
            elseif strcmpi(varargin{i}(1:min([4 end])),'marg')
                indKept = ones(size(listBan,2),1);
                for j = 1:size(listBan,2)
                    if strcmpi(listBan{j},'dataMarg')
                        indKept(j) = 0;
                    end
                    if strcmpi(listBan{j},'dataMargFull')
                        indKept(j) = 0;
                    end
                end
                listBan = listBan(indKept==1);
            else
                warning("\nThe input #%d '%s' was not identified hence ignored. Please try again if necessary.\n",i,varargin{i})
            end
        end
    end

    listVar = evalin("base","who");

    indKept = ones(size(listVar,1),1);
    for i = 1:size(indKept,1)
        for j = 1:size(listBan,2)
            if strcmpi(listBan{j},listVar{i})==1
                indKept(i) = 0;
            end
        end
    end

    listSave = listVar(indKept==1);

    if isempty(listSave)==0
        evalin("base",sprintf("save('%s','%s')",filename,listSave{1}))
        for i = 2:size(listSave,1)
            evalin("base",sprintf("save('%s','%s','-append')",filename,listSave{i}))
        end
    else
        error("\nThere are no variables to save.\n")
    end

end

end

