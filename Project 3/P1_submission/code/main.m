video2images('..\input\detectbuoy.avi',[],'..\input\Images\')

cropImages('..\input\Images\TrainingSet\Frames\','..\input\Images\TrainingSet\CroppedBuoys\')

averageHistogram('..\input\Images\TrainingSet\Frames\','..\input\Images\TrainingSet\CroppedBuoys\','RGB')



