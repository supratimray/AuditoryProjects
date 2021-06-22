function [goodDigitalEvents,goodDigitalTimeStamps] = getGoodDigitalCodes(digitalEvents,digitalTimeStamps,goodEventNumbers)

if ~exist('goodEventNumbers','var'); goodEventNumbers = [9 18 21:52];   end

goodDigitalEvents = [];
goodDigitalTimeStamps = [];

for i=1:length(goodEventNumbers)
    goodPos = find(digitalEvents==goodEventNumbers(i));

    goodDigitalEvents = cat(2,goodDigitalEvents,digitalEvents(goodPos));
    goodDigitalTimeStamps = cat(2,goodDigitalTimeStamps,digitalTimeStamps(goodPos));
end
end