#!/usr/bin/env pwsh

$global:dns_list = aws ec2 describe-instances `
--instance-ids i-0a21f1098dcb1575a `
--query 'Reservations[].Instances[].PublicDnsName'
$global:DNS = $dns_list.Split("")[5]

$global:status_list = aws ec2 describe-instances `
--instance-ids i-0a21f1098dcb1575a `
--query  'Reservations[].Instances[].State.Name'
$global:STATUS = $status_list.Split("")[5]

function Wait-Running {
    while($global:STATUS -notmatch "running")
    {
        Write-Host "Checking if status is running..."

       $wr_status = aws ec2 describe-instances `
        --instance-ids i-0a21f1098dcb1575a `
        --query  'Reservations[].Instances[].State.Name'
        Set-Variable -Name status_list -Value $wr_status -Scope Global
        Set-Variable -Name STATUS -Value $status_list.Split("")[5] -Scope Global

        Write-Host "Status is "$global:STATUS"..."
        Start-Sleep -Seconds 10
    }
    Start-Sleep -Seconds 10

    $wr_dns = aws ec2 describe-instances `
    --instance-ids i-0a21f1098dcb1575a `
    --query 'Reservations[].Instances[].PublicDnsName'
    Set-Variable -Name dns_list -Value $wr_dns -Scope Global
    Set-Variable -Name DNS -Value $dns_list.Split("")[5] -Scope Global
}

function Wait-Stopped {
    while($global:STATUS -notmatch "stopped")
    {
        Write-Host "Checking if status is stopped..."

        $ws_status = aws ec2 describe-instances `
        --instance-ids i-0a21f1098dcb1575a `
        --query  'Reservations[].Instances[].State.Name'
        Set-Variable -Name status_list -Value $ws_status -Scope Global
        Set-Variable -Name STATUS -Value $status_list.Split("")[5] -Scope Global

        Write-Host "Status is "$global:STATUS"..."
        Start-Sleep -Seconds 10
    }
}

if ( $global:STATUS -match "stopped" ) {
    aws ec2 start-instances --instance-ids i-0a21f1098dcb1575a
    Write-Host "Starting..."
    Wait-Running
    ssh -i "C:\Users\hakmo\Documents\j.bravo.pem" ubuntu@$global:DNS
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a
} elseif ( $global:STATUS -match "stopping" ) {
    Write-Host "Waiting for instance to be ready to start..."
    Wait-Stopped
    aws ec2 start-instances --instance-ids i-0a21f1098dcb1575a
    Write-Host "Starting..."
    Wait-Running
    ssh -i "C:\Users\hakmo\Documents\j.bravo.pem" ubuntu@$global:DNS
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a
} else {
    Write-Host "Already on"
    ssh -i "C:\Users\hakmo\Documents\j.bravo.pem" ubuntu@$global:DNS
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a
}