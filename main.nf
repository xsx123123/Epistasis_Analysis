#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------- //
//    Workflow
// ---------------- //

workflow {
    // Analysis init input channel
    ch_samples = Channel.fromPath(params.sample_csv)
    .splitCsv(header:true)
    .map { row -> 
        tuple(
            row.sample_id, 
            row.group + '_1',
            file("${params.raw_data_path}/${row.fq1}"), 
            file("${params.raw_data_path}/${row.fq2}")
        )
    }

    ch_samples.view()

    // 错误！sendMail 不能在这里调用
    // 您必须把它移到下面的 onComplete 事件中
}

// ---------------- //
//  Event Handlers
// ---------------- //

// ✅ 正确的做法：在工作流成功完成时发送邮件
workflow.onComplete {
    
    // 'workflow' 变量在这里是可用的
    sendMail {
        to 'you@gmail.com'
        from 'me@gmail.com'
        attach '/some/path/attachment/file.txt'
        attach '/other/path/image.png'
        
        // 您可以动态地生成邮件主题
        subject "✅ Workflow complete: ${workflow.runName}"

        // 邮件正文
        '''
        Hi there,
        
        The Nextflow workflow '${workflow.runName}' completed successfully.
        
        Look! Multi-lines
        mail content!
        ''' 
    }
}

// 💡 最佳实践：也为失败情况添加一个处理器
workflow.onError {
    
    sendMail {
        to 'you@gmail.com'
        from 'me@gmail.com'
        subject "❌ Workflow FAILED: ${workflow.runName}"

        '''
        Hi there,
        
        The workflow '${workflow.runName}' failed with an error.
        
        Error report: ${workflow.errorReport}
        '''
    }
}