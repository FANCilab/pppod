function info = parseScanImageHeader(h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%ensure inputs or each element of input array is a char
if iscell(h)
  h = fun.map(@char, h);
else
  h = char(h);
end

%parse state.level1.level2.level...=value<newline>
%classes are:
%empties, e.g. 'dottednames=[]'
%string values, e.g. 'dottednames='some value'
%numerical values, e.g. 'dottednames=-13.3'
namesOfEmpties = regexp(h, 'state.(?<name>(\w|\.)+)\=\[\][\r\n]', 'names');
namesToStrings = regexp(h, 'state.(?<name>(\w|\.)+)\=''(?<value>.*?)''[\r\n]', 'names');
namesToNums = regexp(h, 'state.(?<name>(\w|\.)+)\=(?<value>\-?\d.*?|\-?Inf|NaN)[\r\n]', 'names');

%parse all
% namesToVals = regexp(h, 'state.(?<name>(\w|\.)+)\=(?<value>.*?)[\r\n]', 'names')

  function s = toStruct(empties, strings, nums)
    s = struct;
    s = structAssign(s, {empties.name}, []);
    s = structAssign(s, {strings.name}, {strings.value});
    values = {nums.value};
    values = sscanf(sprintf('%s#', values{:}), '%g#')';
    s = structAssign(s, {nums.name}, num2cell(values));
  end

if iscell(h)
  info = cellfun(@toStruct, namesOfEmpties, namesToStrings, namesToNums);
else
  info = toStruct(namesOfEmpties, namesToStrings, namesToNums);
end

end

