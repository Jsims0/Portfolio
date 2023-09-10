import re

from base_classes.BaseModules import BSSHAppModule
from utilities.exceptions import ModuleParameterError
from utilities.bssh_utilities import get_project_id, get_biosample_id
from base_classes.XXXXXXXXXXXX import BSSHFile


class NormaliseRegionsModule(BSSHAppModule):
    def _build_app_config(self):
        workgroup = self._get_parameter('workgroup')

        try:
            output_project_id = self._get_parameter('project_id')
        except ModuleParameterError:
            output_project_id = get_project_id(
                self._get_parameter('project_name'), workgroup)

        self.input_name = re.sub(" ", "_", self.task_id)


        self.app_config = {"app_name": self._get_parameter('app_name'),
                           "app_version": self._get_parameter('app_version'),
                           "app_label": self.app_label,
                           "workgroup": self._get_parameter('workgroup')}

        bam = self._collect_input('bam')
        self.app_config['app_parameters'] = {"project-id": output_project_id,
                                             "app-session-name": self.app_label,
                                             "output-bam": self.input_name + '_norm.bam'}

        try:
            project_id = self._get_parameter('project_id')
        except ModuleParameterError:
            project_id = get_project_id(self._get_parameter('project_name'),
                                        self._get_parameter('workgroup'))


    def _set_output(self):
        biosample_name = self._get_parameter('biosample_name')
        target_x = self._get_parameter('target-x')

        self._add_local_file_output("dataset_content",
                                    "output_dataset_content.json")
        self._add_local_file_output("dataset_info", "output_dataset_info.json")
        self._add_local_file_output("dataset_attributes",
                                    "output_dataset_attributes.json")

        self._add_bssh_file_output("bam", "{}_norm.bam$".format(self.input_name))
        self._add_bssh_file_output("bai", "{}.bai$".format(self.input_name))
        self._add_bssh_file_output("tsv", "{}.tsv$".format(self.input_name))
        self._add_bssh_file_output("Makefile", "Makefile")
