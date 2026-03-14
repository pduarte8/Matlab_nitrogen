
%OutputFile = 'Barents_and_Arctic_European_Sector';
%StructurName = 'BarentsSector'
%FileName = 'data_from_GLODAPv2.2023_European_sector_70-80N.nc';
%StructurName = 'EuropeanSector'
%FileName = 'data_from_GLODAPv2.2023_European_sector_80-90N.nc';

%OutputFile = 'Siberian_Sector';
%StructurName = 'SiberianSector7080N'
%FileName = 'data_from_GLODAPv2.2023_Siberian_sector_70-80N.nc';
%StructurName = 'SiberianSector8090N'
%FileName = 'data_from_GLODAPv2.2023_Siberian_sector_80-90N.nc';

%OutputFile = 'Pacific_Sector_E';
%StructurName = 'PacificSector7080N'
%FileName = 'data_from_GLODAPv2.2023_Pacific_sector_70-80N_170-180E.nc';
%StructurName = 'PacificSector8090N'
%FileName = 'data_from_GLODAPv2.2023_Pacific_sector_80-90N_170-180E.nc';

%OutputFile = 'Pacific_Sector_W';
%StructurName = 'PacificSector7080N'
%FileName = 'data_from_GLODAPv2.2023_Pacific_sector_70-80N_150-180W.nc';
%StructurName = 'PacificSector8090N'
%FileName = 'data_from_GLODAPv2.2023_Pacific_sector_80-90N_150-180W.nc';

OutputFile = 'Canadian_Sector';
%StructurName = 'CanadianSector7080N'
%FileName = 'data_from_GLODAPv2.2023_Canadian_sector_70-80N.nc'
StructurName = 'CanadianSector8090N'
FileName = 'data_from_GLODAPv2.2023_Canadian_sector_80-90N.nc'

Vars = {'PRESSURE', 'DEPTH', 'TEMPERATURE','SALNTY','OXYGEN','PHSPHT','SILCAT','NITRAT','NITRIT',...
        'ALKALI','TCARBN','TOC','DOC','DON','TDN','CHLORA','THETA','SIGMA0','AOU'}
%Vars = {'NITRAT'}
for i = 1:length(Vars)
    VarName = string(Vars(i));
    VarFlagName = strcat(VarName,'_qc');
    ncid = netcdf.open(FileName,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,VarFlagName);
    flags = netcdf.getVar(ncid,varid,'int8'); % returns raw int8 values
    netcdf.close(ncid);
    info = ncinfo(FileName);
    allName = {info.Variables.Name};
    idx = find(strcmp(allName, VarName));
    units = info.Variables(idx).Attributes(2).Value; 
    info = ncinfo(FileName, VarFlagName);
    flag_type = info.Attributes(strcmp({info.Attributes.Name}, 'comment')).Value
    
   
    if contains(flag_type, 'ODV', 'IgnoreCase', true)
       if any(flags > 40 & flags < 60)
          flags = flags - 48;
       end
       good_mask = ismember(flags,[0 1]);
    elseif contains(flag_type, 'gtspp', 'IgnoreCase', true)
       if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
       end 
       good_mask = ismember(flags,1);  
    elseif contains(flag_type, 'argo')
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);
    elseif contains(flag_type, 'SeaDataNet')
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);
    elseif contains(flag_type, 'ESEAS')
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);  
    elseif contains(flag_type, 'WOD') % includes WODSTATION
        if any(flags > 40 & flags < 60)
          flags = flags - 48;
        end 
        good_mask = ismember(flags, 0); 
    elseif contains(flag_type, 'WOCE', 'IgnoreCase', true) % includes WOCEBOTTLE, WOCECTD, WOCESAMPLE
        if any(flags > 40 & flags < 60)
          flags = flags - 48;
        end  
        good_mask = ismember(flags, 2);
    elseif contains(flag_type, 'QARTOD', 'IgnoreCase', true)
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);
    elseif contains(flag_type, 'OceanSITES', 'IgnoreCase', true)
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);
    elseif contains(flag_type, 'IODE', 'IgnoreCase', true)
        if any(flags > 40 & flags < 60)
          flags = flags + 1 - 48;
        end 
        good_mask = ismember(flags, 1);    
    else
        good_mask = []; % fallback
    end
    if i == 1
    info = ncinfo(FileName);
    S.(StructurName).longitude.Value = ncread(FileName, 'longitude');
    idx = find(strcmp(allName, 'longitude'));
    units = info.Variables(idx).Attributes(3).Value;  
    S.(StructurName).longitude.Units = units; 
    S.(StructurName).latitude.Value = ncread(FileName, 'latitude');
    idx = find(strcmp(allName, 'latitude'));
    units = info.Variables(idx).Attributes(3).Value; 
    S.(StructurName).latitude.Units = units;
    S.(StructurName).date_time.Value = ncread(FileName, 'date_time');
    idx = find(strcmp(allName, 'date_time'));
    units = info.Variables(idx).Attributes(3).Value; 
    S.(StructurName).date_time.Units = units;
    seconds = ncread(FileName, 'date_time') * 24 * 3600;
    if strcmp(StructurName,'EuropeanSector7080N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1972);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);
    elseif strcmp(StructurName,'EuropeanSector8090N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1980);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);
    elseif strcmp(StructurName,'SiberianSector7080N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1993);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);   
    elseif strcmp(StructurName,'SiberianSector8090N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1991);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D); 
    elseif strcmp(StructurName,'PacificSector7080N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1993);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);   
    elseif strcmp(StructurName,'PacificSector8090N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1994);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);  
    elseif strcmp(StructurName,'CanadianSector7080N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1974);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);   
    elseif strcmp(StructurName,'CanadianSector8090N')
       [Y M D] = unixsecs2date_from_any_initial_year(seconds,1996);
       UNIXseconds = date2unixsecs_from_any_initial_year(1970,Y,M,D);     
    end
    S.(StructurName).UNIXSECS.Value = UNIXseconds;
    S.(StructurName).UNIXSECS.Units = 'Unix seconds since 1970-01-01 00:00:00 UTC';
    S.(StructurName).Year.Value = Y;
    S.(StructurName).Month.Value = M;
    S.(StructurName).Day.Value = D;
    end
    Data = ncread(FileName, VarName);
    Data(~good_mask) = NaN;
    S.(StructurName).(VarName).Value = Data;
    S.(StructurName).(VarName).Units = units;
end
isfile([OutputFile,'.mat'])
if ~isfile([OutputFile,'.mat'])
    save(OutputFile,'S')
else    
    save(OutputFile,'S','-append')
end