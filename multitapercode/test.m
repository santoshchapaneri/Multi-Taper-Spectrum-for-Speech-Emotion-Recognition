frameNum = 60;

subplot(411);plot(spec_hamming(:,frameNum));
subplot(412);plot(spec_thomson(:,frameNum));
subplot(413);plot(spec_multip(:,frameNum));
subplot(414);plot(spec_SWCE(:,frameNum));