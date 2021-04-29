function testView(){
    $("#log_in_box").hide();
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
    $("#test-env").show();
}

function switchToLogInView(){
    $("#userNameInput").val('');
    $("#passwordInput").val('');
    $("#log_in_box").show();
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
    $("#email").hide();
}

function switchToLoggedInView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
}

function switchToPipelineView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
    document.getElementById("runs_list").innerHTML = '';
    document.getElementById("run_button").remove();
}

function switchToTypeView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
    document.getElementById("runs_list").innerHTML = '';
    document.getElementById("run_button").remove();
}

function switchToRunIdView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
}

function switchToLaneView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
}

function switchToMatchView(){
    $("#log_in_box").hide();
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
    $("#email").hide();
}

function switchToReportView(){
    $("#log_in_box").hide();
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
    document.getElementById("report_runs_list").innerHTML = '';
    document.getElementById("samples").innerHTML = '';
    // document.getElementById("panels").innerHTML = '';
    get_report_runs();
}

function switchToDocumentationView(){
    window.location.replace("toc.html")
}