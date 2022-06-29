function task = HCPtaskdef(mat_files)
task=[];
task.Wrkmem.TIM.defining_event = 'onset of an image that the subject has to match or not with the target image';
task.Wrkmem.TIM.trial_duration = [-1.5 2.5];
task.Wrkmem.TIM.studied = {'networks employed during image processing (faces vs tools)', 'networks of memory recall for different memory loads'};
task.Wrkmem.TIM.index = (contains(mat_files, 'Wrkmem') + contains(mat_files, 'TIM')) == 2;
task.Wrkmem.TIM.cycletime = {};
task.Wrkmem.TIM.lifetimes = {};

task.Wrkmem.TRESP.defining_event = '';
task.Wrkmem.TRESP.trial_duration = [-1.5 1.5];
task.Wrkmem.TRESP.studied = {'mechanism of memory recall and decision making just prior to the response mechanism of memory storage after response when the memory load is dynamic (2-back) as compared to static', 'motor network'};
task.Wrkmem.TRESP.index = (contains(mat_files, 'Wrkmem') + contains(mat_files, 'TRESP')) == 2;
task.Wrkmem.TRESP.cycletime = {};
task.Wrkmem.TRESP.lifetimes = {};

task.StoryM.TEV.defining_event = 'onset of any task event during the question and option period in stories and math problems';
task.StoryM.TEV.trial_duration = [-1.5 1.5];
task.StoryM.TEV.studied = {'networks employed during math calculations', 'compare networks during first parts of a calculation and later parts with heavier load', 'networks during the onset of a story sentence relative to that of a math problem sentence', 'networks during the presentation of the correct answer as compared to the wrong one', 'activated areas when a number is presented as compared to when an operand is presented'};
task.StoryM.TEV.index = (contains(mat_files, 'StoryM') + contains(mat_files, 'TEV')) == 2;
task.StoryM.TEV.cycletime = {};
task.StoryM.TEV.lifetimes = {};


task.StoryM.TRESP.defining_event = 'onset of button press by the subject';
task.StoryM.TRESP.trial_duration = [-1.5 1.5];
task.StoryM.TRESP.studied = {'motor network when a finger has already been picked for response', 'compare immediate responses with later responses which can mean uncertainty in the answer'};
task.StoryM.TRESP.index = (contains(mat_files, 'StoryM') + contains(mat_files, 'TRESP')) == 2;
task.StoryM.TRESP.cycletime = {};
task.StoryM.TRESP.lifetimes = {};

task.StoryM.BSENT.defining_event = 'trials contain entire sentences. For stories, this is a sentence during narration without the option sentence at the end of the story. For math, this is the sentence of the math problem excluding the option sentence at the end.';
task.StoryM.BSENT.trial_duration = [-1 inf]; % till 1 sec after its offset
task.StoryM.BSENT.studied = {'compare how brain networks change at the later parts of the story as compared to earlier parts of a story', 'compare how brain networks change at the later math problems in a block as compared to earlier math problems in a block'};
task.StoryM.BSENT.index = (contains(mat_files, 'StoryM') + contains(mat_files, 'BSENT')) == 2;
task.StoryM.BSENT.cycletime = {};
task.StoryM.BSENT.lifetimes = {};

task.StoryM.BUN.defining_event = 'rials containing entire Blocks of stimulus Units As stimulus unit is defined an entire story or an entire math problem including the option part.';
task.StoryM.BUN.trial_duration = [-1 inf]; % till 1 sec after its offset
task.StoryM.BUN.studied = {'compare the overall networks deployed during the processing of stories vs the processing of math problems', 'The logic behind selecting entire units as trials comes from the fact that the total duration of stories is balanced with the total duration of math problems. However the number of math problems is significantly higher than the number of stories and the inverse holds for their duration. So actual comparison between stories and math can be performed at the block level and not at event level.'};
task.StoryM.BUN.index = (contains(mat_files, 'StoryM') + contains(mat_files, 'BUN')) == 2;
task.StoryM.BUN.cycletime = {};
task.StoryM.BUN.lifetimes = {};