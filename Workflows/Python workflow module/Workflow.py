from base_classes.XXXXXXXXXXX import XXXXXXXXXXXXXXX
from modules.MonitorRunModule import MonitorRunModule
from modules.PadHopperModule import PadHopperModule
from modules.DragenGermlineModule import DragenGermlineModule
from modules.NormaliseRegionsModule import NormaliseRegionsModule
from modules.XXXXXXXXXXXXX import XXXXXXXXXXXXXX

class XXXXXXXWorkflow(Workflow):

    def setup_workflow(self):
        bs_set_name = [*self.analysis_config['runs'].keys()][0]

        bs = self.analysis_config['runs'][bs_set_name]
        ## DRAGEN alignment task for each biosample
        alignment_task = DragenGermlineModule(
            module_id='dragen_alignment',
            task_id=bs_set_name,
            task_parameters={'biosample_names': bs},
            workflow_config=self.workflow_config,
            out_dir=self.out_dir,
            existing_app_sessions=self.existing_app_sessions)

        #alignment_task.add_dependency(monitor_run_task)
        self._add_task(alignment_task)

        ## normailse regions task
        nr_task = NormaliseRegionsModule(
            module_id='normalize_regions',
            task_id="{}X_{}".format(
                self._get_parameter('target-x'), bs_set_name),
            task_parameters={'biosample_name': bs },
            workflow_config=self.workflow_config,
            out_dir=self.out_dir,
            existing_app_sessions=self.existing_app_sessions)

        nr_task.add_dependency(alignment_task)
        self._add_task(nr_task)

        ## Full DRAGEN germline analysis task
        germline_task = DragenGermlineModule(
            module_id='dragen_germline',
            task_id="{}X_{}".format(
                self._get_parameter('target-x'), bs_set_name),
            task_parameters={'biosample_name': bs},
            workflow_config=self.workflow_config,
            out_dir=self.out_dir,
            existing_app_sessions=self.existing_app_sessions)

        germline_task.add_dependency(nr_task)
        self._add_task(germline_task)

        ## XXXXXXXXX analysis task
        XXXXXXXXXX_task = XXXXXXXXXXModule(
            module_id='XXXXXXXXXX',
            task_id="{}_{}X".format(
                self._get_parameter('target-x'), bs_set_name),
            task_parameters={'biosample_name': bs},
            workflow_config=self.workflow_config,
            out_dir=self.out_dir,
            existing_app_sessions=self.existing_app_sessions)
        firebrand_task.add_dependency(germline_task)
        self._add_task(XXXXXXXX_task)
