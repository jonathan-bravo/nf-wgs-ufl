function show_login() {
    $('.log-in-box-shadow').css('opacity', '1');
    $('.log-in-box-shadow').removeClass('inactive');
    $('.log-in-box-shadow').addClass('active');
    $('.log-in-box-shadow').css('zIndex','-1');
    setTimeout(function() {
        $('.log-in-box').css('opacity', '1');
        $('.log-in-box').css('zIndex','1');
        $('.log-in-box').css('box-shadow', '10px 10px 5px rgba(0, 0, 0, 0.1)');
        $('.log-in-box').removeClass('inactive');
        $('.log-in-box').addClass('active');
    }, 300);

    $('.request-access-box').removeClass('active');
    $('.request-access-box').addClass('inactive');
    setTimeout(function(){
        $('.request-access-box-shadow').removeClass('active');
        $('.request-access-box-shadow').addClass('inactive');
    }, 300);
    setTimeout(function(){
        $('.request-access-box-shadow').css('opacity', '0');
        $('.request-access-box-shadow').css('zIndex', '-1');
        $('.request-access-box').css('opacity', '0');
        $('.request-access-box').css('box-shadow', 'none');
        $('.request-access-box').css('zIndex', '-1');
    }, 1500);
}

function show_request_access() {

    $('.request-access-box-shadow').css('opacity', '1');
    $('.request-access-box-shadow').removeClass('inactive');
    $('.request-access-box-shadow').addClass('active');
    $('.request-access-box-shadow').css('zIndex','-1');
    setTimeout(function() {
        $('.request-access-box').css('opacity', '1');
        $('.request-access-box').css('zIndex', '1');
        $('.request-access-box').css('box-shadow', '10px 10px 5px rgba(0, 0, 0, 0.1)');
        $('.request-access-box').removeClass('inactive');
        $('.request-access-box').addClass('active');
    }, 300);

    $('.log-in-box').removeClass('active');
    $('.log-in-box').addClass('inactive');
    setTimeout(function(){
        $('.log-in-box-shadow').removeClass('active');
        $('.log-in-box-shadow').addClass('inactive');
    }, 300);
    setTimeout(function(){
        $('.log-in-box-shadow').css('opacity', '0');
        $('.log-in-box-shadow').css('zIndex','-1');
        $('.log-in-box').css('opacity', '0');
        $('.log-in-box').css('box-shadow', 'none');
        $('.log-in-box').css('zIndex','-1');
    }, 1500);
}

function testView(){
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    $("#request_runs_box").hide();
    $("#test-env").show();
}

function switchToLogInView(){
    $("#userNameInput").val('');
    $("#passwordInput").val('');
    $("#intro").show();
    $("#log_out_button").hide();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
    $("#transfer_run_box").hide();
}

function switchToLoggedInView(){
    $("#intro").hide();
    $("#menu").show();
    $("#log_out_button").show();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
}

function switchToPipelineView(){
    $("#intro").hide();
    $("#logOutButton").show();
    $("#menu").hide();
    $("#workflow_question_box").show();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
    document.getElementById("runs_list").innerHTML = '';
    document.getElementById("run_button").remove();
}

function switchToTypeView(){
    $("#intro").hide();
    $("#logOutButton").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").show();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
    document.getElementById("runs_list").innerHTML = '';
    document.getElementById("run_button").remove();
}

function switchToRunIdView(){
    $("#intro").hide();
    $("#logOutButton").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").show();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
}

function switchToLaneView(){
    $("#intro").hide();
    $("#logOutButton").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").show();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
}

function switchToMatchView(){
    $("#intro").hide();
    $("#logOutButton").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").show();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report_runs_box").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
}

function switchToReportView(){
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#request_runs_box").hide();
    $("#email").hide();
    document.getElementById("report_runs_list").innerHTML = '';
    document.getElementById("samples").innerHTML = '<div id="parent_selector_box" class="flex-child"><div id="parent_selector_child_box" class="flex-child-child"><select multiple="multiple" id="parent_selector"></select></div></div>';
    get_report_runs();
}

function switchToRequestView() {
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    document.getElementById("request_runs_list").innerHTML = '';
    get_request_runs();
}

function switchToRunOrSampleView() {
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    $("#request_runs_box").hide();
    $("#request_run_or_sample_box").show();
}

function switchToRequestSamplesView() {
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    $("#request_runs_box").hide();
    $("#request_run_or_sample_box").hide();
    document.getElementById("request_samples").innerHTML = '';
    $("#request_samples_box").show();
}

function switchToRequestReportsView() {
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    $("#request_runs_box").hide();
    $("#request_run_or_sample_box").hide();
    $("#request_samples_box").hide();
    $("#request_reports_box").show();
}

function switchToTransferView() {
    $("#intro").hide();
    $("#log_out_button").show();
    $("#menu").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#germline_command").hide();
    $("#multiqc_command").hide();
    $("#report_selection_back").hide();
    $("#report_selection").hide();
    $("#email").hide();
    $("#request_runs_box").hide();
    $("#request_run_or_sample_box").hide();
    $("#request_samples_box").hide();
    $("#request_reports_box").hide();
    $("#transfer_run_box").show();
}

function switchToDocumentationView(){
    window.location.replace("toc.html")
}