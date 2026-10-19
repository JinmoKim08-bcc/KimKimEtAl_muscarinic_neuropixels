function orthonormal_basis = Perform_CD_analysis_v4b(CR, CL, ER, EL, param)

RL = [CR CL];
RL = RL - repmat(mean(RL, 2), 1, size(RL, 2)); % remove mean
[u s v] = svd(RL');

wt = (CR-EL)/2 + (ER-CL)/2;
i_t = (param.ITI_time*param.fs_bin)+1:((param.ITI_time+param.sample_time)*param.fs_bin);
CD_stim_mode = mean(wt(:,i_t),2);

i_t = floor(((param.ITI_time+param.sample_time)*param.fs_bin)+1:((param.ITI_time+param.sample_time+0.3)*param.fs_bin));
CD_stim_mode_after = mean(wt(:,i_t),2);

wt = (CR-ER)/2 + (EL-CL)/2;
i_t = ((param.ITI_time+param.sample_time)*param.fs_bin+1):((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_choice_mode = mean(wt(:,i_t),2);
wt = (CR-ER)/2 + (CL-EL)/2;
i_t = ((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1):((param.ITI_time+param.sample_time+param.delay_time+param.response_time)*param.fs_bin);
CD_outcome_mode = mean(wt(:,i_t),2);

wt = (CR-CL)/2;
i_t = floor(((param.ITI_time+0.01)*param.fs_bin)+1):floor(((param.ITI_time+0.05)*param.fs_bin));
CD_sample_mode = mean(wt(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+param.delay_time-0.6)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_delay_mode = mean(wt(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time+0.4)*param.fs_bin);
CD_go_mode = mean(wt(:,i_t),2);

wt = (CR+CL)/2;
i_t1 = ((param.ITI_time-0.6)*param.fs_bin+1):(param.ITI_time*param.fs_bin);
i_t2 = floor((param.ITI_time+param.sample_time+param.delay_time-0.1)*param.fs_bin+1):((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
Ramping_mode = mean(wt(:,i_t2),2)-mean(wt(:,i_t1),2);
i_t1 = floor(((param.ITI_time+param.sample_time+param.delay_time-0.1)*param.fs_bin)+1):(floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin));
i_t2 = floor(((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1)):floor(((param.ITI_time+param.sample_time+param.delay_time+0.1)*param.fs_bin));
GoDirection_mode = mean(wt(:,i_t2),2)-mean(wt(:,i_t1),2);

wt = (CR-CL)/2;
i_t = floor((param.ITI_time+param.sample_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_delay_mode_rah = mean(wt(:,i_t),2);

wt = (CR-CL)/2;
i_t = floor((param.ITI_time+param.sample_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+0.66)*param.fs_bin);
CD_delay_1mode_rah = mean(wt(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+0.66)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+0.66+0.66)*param.fs_bin);
CD_delay_2mode_rah = mean(wt(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+0.66+0.66)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_delay_3mode_rah = mean(wt(:,i_t),2);


wt_jm = ((CR+CL)-(ER+EL))/2;
i_t = floor(((param.ITI_time+0.01)*param.fs_bin)+1):floor(((param.ITI_time+0.05)*param.fs_bin));
CD_sample_mode_jm = mean(wt_jm(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+param.delay_time-0.6)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_delay_mode_jm = mean(wt_jm(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time+0.4)*param.fs_bin);
CD_go_mode_jm = mean(wt_jm(:,i_t),2);

wt_jm = ((CR-ER)+(CL-EL))/2;
i_t1 = ((param.ITI_time-0.6)*param.fs_bin+1):(param.ITI_time*param.fs_bin);
i_t2 = floor((param.ITI_time+param.sample_time+param.delay_time-0.1)*param.fs_bin+1):((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
Ramping_mode_jm = mean(wt_jm(:,i_t2),2)-mean(wt_jm(:,i_t1),2);
i_t1 = floor(((param.ITI_time+param.sample_time+param.delay_time-0.1)*param.fs_bin)+1):(floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin));
i_t2 = floor(((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1)):floor(((param.ITI_time+param.sample_time+param.delay_time+0.1)*param.fs_bin));
GoDirection_mode_jm = mean(wt_jm(:,i_t2),2)-mean(wt_jm(:,i_t1),2);


%%%%%%% simple %%%%%%%
wt = (CR-CL)/2;
i_t = floor(((param.ITI_time)*param.fs_bin)+1):floor(((param.ITI_time+param.sample_time)*param.fs_bin));
CD_sample_mode_simple = mean(wt(:,i_t),2);
wt = (CR-CL)/2;
i_t = floor((param.ITI_time+param.sample_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin);
CD_delay_mode_simple = mean(wt(:,i_t),2);
i_t = floor((param.ITI_time+param.sample_time+param.delay_time)*param.fs_bin+1):floor((param.ITI_time+param.sample_time+param.delay_time+0.5)*param.fs_bin);
CD_response_mode_simple = mean(wt(:,i_t),2);





%%%%%%%
CD_stim_mode = CD_stim_mode / norm(CD_stim_mode);
CD_choice_mode = CD_choice_mode / norm(CD_choice_mode);
CD_outcome_mode = CD_outcome_mode / norm(CD_outcome_mode);
CD_sample_mode = CD_sample_mode / norm(CD_sample_mode);
CD_delay_mode = CD_delay_mode / norm(CD_delay_mode);
CD_delay_mode_rah = CD_delay_mode_rah / norm(CD_delay_mode_rah);
CD_delay_1mode_rah = CD_delay_1mode_rah / norm(CD_delay_1mode_rah);
CD_delay_2mode_rah = CD_delay_2mode_rah / norm(CD_delay_2mode_rah);
CD_delay_3mode_rah = CD_delay_3mode_rah / norm(CD_delay_3mode_rah);
CD_go_mode = CD_go_mode / norm(CD_go_mode);
Ramping_mode = Ramping_mode / norm(Ramping_mode);
GoDirection_mode = GoDirection_mode / norm(GoDirection_mode);
CD_sample_mode_jm = CD_sample_mode_jm / norm(CD_sample_mode_jm);
CD_delay_mode_jm = CD_delay_mode_jm / norm(CD_delay_mode_jm);
CD_go_mode_jm = CD_go_mode_jm / norm(CD_go_mode_jm);
Ramping_mode_jm = Ramping_mode_jm / norm(Ramping_mode_jm);
GoDirection_mode_jm = GoDirection_mode_jm / norm(GoDirection_mode_jm);
CD_stim_mode_after = CD_stim_mode_after / norm(CD_stim_mode_after);

CD_sample_mode_simple = CD_sample_mode_simple / norm(CD_sample_mode_simple);
CD_delay_mode_simple = CD_delay_mode_simple / norm(CD_delay_mode_simple);
CD_response_mode_simple = CD_response_mode_simple / norm(CD_response_mode_simple);

if strcmp(param.cd_opt, 'yang')
    orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_mode CD_go_mode Ramping_mode GoDirection_mode v]);
elseif strcmp(param.cd_opt, 'inagaki')
    orthonormal_basis = Gram_Schmidt_process([CD_sample_mode CD_delay_mode CD_go_mode Ramping_mode GoDirection_mode v]);
%     orthonormal_basis = Gram_Schmidt_process([CD_delay_mode CD_go_mode Ramping_mode GoDirection_mode v]);
elseif strcmp(param.cd_opt, 'jinmo')
%     orthonormal_basis = Gram_Schmidt_process([CD_sample_mode_jm CD_delay_mode_jm CD_go_mode_jm Ramping_mode_jm GoDirection_mode_jm v]);
%     orthonormal_basis = Gram_Schmidt_process([CD_sample_mode_jm CD_delay_mode_jm CD_go_mode_jm v]);
    orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_mode CD_go_mode Ramping_mode GoDirection_mode v]);
elseif strcmp(param.cd_opt, 'rah')
    orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_mode_rah CD_go_mode Ramping_mode GoDirection_mode v]);
elseif strcmp(param.cd_opt, 'rah_3delay')
    orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_1mode_rah CD_delay_2mode_rah CD_delay_3mode_rah CD_go_mode Ramping_mode GoDirection_mode v]);
elseif strcmp(param.cd_opt, 'test')
    orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_mode CD_go_mode Ramping_mode GoDirection_mode CD_stim_mode_after v]);
elseif strcmp(param.cd_opt, 'simple')
    orthonormal_basis = Gram_Schmidt_process([CD_sample_mode_simple CD_delay_mode_simple CD_response_mode_simple v]);
elseif strcmp(param.cd_opt, 'simple_w_ramping_go')
    orthonormal_basis = Gram_Schmidt_process([CD_sample_mode_simple CD_delay_mode_simple CD_response_mode_simple Ramping_mode GoDirection_mode v]);
end

end

