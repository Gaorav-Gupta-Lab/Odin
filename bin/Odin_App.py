#!/usr/bin/env python3
"""

Code for a GUI to generate Odin script file.
@author: Dennis A. Simpson
         RTP Genomics LLC
         Chapel Hill, NC  27599
@copyright: 2019
"""
import datetime
import os
import collections
import sys
from contextlib import suppress
from pathlib import Path
import subprocess
import dill
import wx
import wx.adv
from wx.lib.wordwrap import wordwrap
import wx.lib.sized_controls as sized_controls
import wx.lib.intctrl

__author__ = 'Dennis A. Simpson'
__version__ = '0.4.0'
__package__ = 'Odin'
__copyright__ = '(C) 2019'


class IntBoxes:
    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_int_caller(self, name, value=None, validator=None):
        """

        :param validator:
        :param name:
        :param value:
        :return:
        """

        int_ctrl = wx.lib.intctrl.IntCtrl(self.main_panel, wx.ID_ANY, name=name, value=value, allow_none=1, min=0)
        if validator == "Number":
            int_ctrl.SetValidator(NumValidator(self.main_panel, name))
        int_ctrl.SetFont(wx.Font(11, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))
        int_ctrl.SetColors(default_color=wx.BLACK, oob_color=wx.RED)
        return int_ctrl


class DataBoxes:
    """
    Generic ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_box(self, data_list, name, value=None):
        """

        :param data_list:
        :param name:
        :param value:
        :return:
        """
        if value is None:
            value = ""
        my_combo = wx.ComboBox(self.main_panel, wx.ID_DEFAULT, value, choices=data_list, style=wx.CB_SORT, name=name)
        my_combo.SetFont(wx.Font(11, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))
        return my_combo


class MaskedDataBoxes:
    """
    Generic Masked ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_masked_box(self, data_list, name, default_value=None, validator=None):
        """

        :param validator:
        :param data_list:
        :param name:
        :param default_value:
        :return:
        """

        my_masked_combo = \
            wx.ComboBox(self.main_panel, id=-1, choices=data_list, style=wx.CB_SORT, name=name)

        if validator == "Boolean":
            my_masked_combo.SetValidator(BooleanValidator(self.main_panel, name))

        if default_value:
            my_masked_combo.SetValue(default_value)
        my_masked_combo.SetFont(wx.Font(11, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))

        return my_masked_combo


class ButtonGenerator:
    """
    Generic Button Generator for Project Insuring all Buttons are the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_buttons(self, label):
        my_button = wx.Button(self.main_panel, wx.ID_ANY, label=label)
        my_button.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))
        return my_button


class BooleanValidator(wx.Validator):
    def __init__(self, parent, name):
        super(BooleanValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        return BooleanValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if text in ["True", "False", "Mouse", "Human"]:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            if self.name == "--Species":
                message = "{} allowed values are 'Mouse' or 'Human'".format(self.name)
            else:
                message = "{} values are 'True' or 'False'".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class PathValidator(wx.Validator):
    def __init__(self, parent, name):
        super(PathValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return PathValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if "/" in text:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} requires full path statements".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class TextValidator(wx.Validator):
    def __init__(self, parent, name):
        super(TextValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return TextValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if len(text) > 0:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} cannot be Null".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
        return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class NumValidator(wx.lib.intctrl.IntValidator):
    def __init__(self, parent, name):
        wx.lib.intctrl.IntValidator.__init__(self)
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return NumValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()
        if isinstance(text, int):
            if (self.name == "--Spawn" and text > 0) or self.name != "--Spawn" and text >= 0:
                textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                textCtrl.Refresh()
                return True
        else:
            if self.name == "--Spawn":
                message = "{} takes integers >0 only".format(self.name)
            else:
                message = "{} takes positive integers only".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()

            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class CommonControls:
    def __init__(self, parent):
        self.parent = parent
        self.default_data = wx.GetTopLevelParent(parent).default_data
        self.dataframe = wx.GetTopLevelParent(parent).dataframe
        self.button_generator = ButtonGenerator(parent)
        self.databoxes = DataBoxes(parent)
        self.my_masked_databoxes = MaskedDataBoxes(parent)
        self.my_int_box = IntBoxes(parent)

    def add_line(self):
        label = wx.StaticText(self.parent, wx.ID_DEFAULT, "", style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))

        label.SetMinSize((200, 10))
        line = wx.StaticLine(self.parent, wx.ID_ANY, style=wx.LI_HORIZONTAL, name="line")
        line.SetBackgroundColour('#E5E8E8')

        return label, line

    def folder_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        working_dir_btn = self.button_generator.my_buttons(name)
        working_dir_btn.Bind(wx.EVT_BUTTON, lambda event, temp=working_dir_btn: self.select_folder(event, input_ctrl))

        return working_dir_btn, input_ctrl

    def file_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        fastq_file_btn = self.button_generator.my_buttons(name)
        fastq_file_btn.Bind(wx.EVT_BUTTON, lambda event, temp=fastq_file_btn: self.select_file(event, input_ctrl))

        return fastq_file_btn, input_ctrl

    def restricted_selector(self, name, value=None, validator=None, tip=None):
        """
        Restricts input.
        :param validator:
        :param name:
        :param value:
        :param tip:
        :return:
        """
        input_ctrl = self.my_masked_databoxes.my_masked_box(self.default_data[name], name, value, validator)
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        label = wx.StaticText(self.parent, wx.ID_ANY, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def text_control(self, name, tip=None):
        try:
            value = self.default_data[name]
        except KeyError:
            value = ""
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name, value)
        input_ctrl.SetValidator(TextValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        # label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def int_control(self, name, validator=None, tip=None):
        """

        :param validator:
        :param name:
        :param tip:
        :return:
        """

        input_ctrl = self.my_int_box.my_int_caller(name, self.default_data[name], validator)
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        # label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))

        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def select_folder(self, event, ctrl):
        dlg = wx.DirDialog(self.parent, "Choose a Folder")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def select_file(self, event, ctrl):
        dlg = wx.FileDialog(self.parent, "Choose a file")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def int_control_builder(self, name, default_value):
        name = name
        try:
            value = default_value[name]
        except KeyError:
            value = None

        input_ctrl = self.my_int_box.my_int_caller(name, value)

        return input_ctrl


class SkadiErrorMatrixPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SkadiErrorMatrixPanel, self).__init__(parent, wx.ID_ANY, name="SkadiErrorMatrixPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Build Error Matrix")
        # title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Unfiltered VCF File"),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name")]

        # self.panel_grid_sizer.SetEmptyCellSize((10, 5))
        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class SkadiFilterVCFPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SkadiFilterVCFPanel, self).__init__(parent, wx.ID_ANY, name="SkadiFilterVCFPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Skadi Filter VCF")
        # title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.file_selector("--Tumor_VCF", tip="Full Path to Unfiltered VCF File"),
             my_controls.file_selector("--Ref_Seq", tip="Full Path to Genomic Reference Sequence"),
             my_controls.file_selector("--dbSNP", tip="Full Path to dbSNP file"),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name"),
             my_controls.restricted_selector("--Include_All", validator="Boolean",
                                             tip="Include normal positions in VCF file?"),
             my_controls.text_control("--SNP_Error_Freq", tip="For Filtering.  SNP Frequency in Error Matrix"),
             my_controls.text_control("--Minimum_Allele_Freq",
                                      tip="When Filtering; What is the Minimum Minimum ALT Frequency Allowed"),
             my_controls.int_control("--Min_Fold_Increase", validator="Number",
                                     tip="When Filtering; What is the Minimum Increase of the Alt Allele"),
             my_controls.int_control("--Minimum_Allele_Count", validator="Number", tip="Minimum Alt Allele Count"),
             my_controls.int_control("--Minimum_Read_Depth", validator="Number"),
             ]

        # self.panel_grid_sizer.SetEmptyCellSize((10, 5))
        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class AnnotateVCFPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(AnnotateVCFPanel, self).__init__(parent, wx.ID_ANY, name="AnnotateVCFPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Annotate VCF")
        # title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.folder_selector("--snpEff_Folder", tip="Full Path to snpEff install folder"),
             my_controls.file_selector("--Input_VCF", tip="Full Path to Input_VCF file"),
             my_controls.int_control("--Java_Mem", validator="Number", tip="How many Gb memory to give java?"),
             my_controls.text_control("--Dataset", tip="snpEff Dataset to use"),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name"),
             ]

        # self.panel_grid_sizer.SetEmptyCellSize((10, 5))
        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class SkadiPanelVariantSearchPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(SkadiPanelVariantSearchPanel, self).__init__(parent, wx.ID_ANY, name="SkadiPanelVariantSearchPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Skadi Variant Search")
        # title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.file_selector("--BAM_File", tip="Full Path to Sorted and Indexed BAM File"),
             my_controls.file_selector("--Ref_Seq", tip="Full Path to Genomic Reference Sequence"),
             my_controls.file_selector("--dbSNP", tip="Full Path to dbSNP file"),
             my_controls.file_selector("--Target_File", tip="Full Path to Genome Target File"),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.add_line(),

             my_controls.restricted_selector("--Verbose", value="INFO"),
             my_controls.text_control("--Job_Name"),
             my_controls.int_control("--Spawn", validator="Number", tip="How Many Threads/Cores to Run."),
             my_controls.restricted_selector("--Compression_Level", "9"),
             my_controls.restricted_selector("--Species", validator="Boolean"),
             my_controls.restricted_selector("--Demultiplex", validator="Boolean"),
             my_controls.restricted_selector("--Include_All", validator="Boolean", value=False,
                                             tip="Include normal positions in VCF file?"),
             my_controls.text_control("--N_Limit", tip="Fraction of Sequence Allowed to be N"),
             my_controls.restricted_selector("--Target"),
             my_controls.restricted_selector("--Filter_VCF_File", validator="Boolean", value=False),
             my_controls.int_control("--Boundary_Padding", validator="Number",
                                     tip="Number of Additional Nucleotides Outside of Regions"),
             my_controls.restricted_selector("--Strict_Boundaries", validator="Boolean", value="True"),
             my_controls.int_control("--Minimum_Allele_Count", validator="Number", tip="Minimum Alt Allele Count"),
             my_controls.int_control("--Minimum_Read_Depth", validator="Number"),
             my_controls.int_control("--Minimum_Base_Quality", validator="Number")]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class ThruPlexPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(ThruPlexPanel, self).__init__(parent, wx.ID_ANY, name="ThruPlexPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        self.panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        # title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        # title_icon = wx.StaticBitmap(self, wx.ID_DEFAULT, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Odin ThruPLEX Pipeline")
        title.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.FONTWEIGHT_BOLD, faceName='Inconsolata'))
        title.SetMinSize((int(parent.GetSize()[0] * 0.5), 45))
        title.SetForegroundColour("blue")
        title_sizer.Add(title, 0, wx.ALL)
        panel_sizer.Add(title_sizer, 0, wx.ALL | wx.CENTER)
        panel_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.file_selector("--FASTQ1", tip="Full Path to FASTQ1"),
             my_controls.file_selector("--FASTQ2", tip="Full Path to FASTQ2"),
             my_controls.file_selector("--Ref_Seq", tip="Full Path to Genomic Reference Sequence"),
             my_controls.file_selector("--Aligner_Ref_Seq", tip="Full Path to Genomic Aligner Reference Sequence"),
             my_controls.file_selector("--Index_File"),
             my_controls.file_selector("--dbSNP", tip="Full Path to dbSNP file"),
             my_controls.file_selector("--Target_File", tip="Full Path to Genome Target File"),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Verbose", value="INFO"),
             my_controls.text_control("--Job_Name"),
             my_controls.int_control("--Spawn", validator="Number"),
             my_controls.restricted_selector("--Compression_Level", "9"),
             my_controls.restricted_selector("--Species", validator="Boolean",
                                             tip="Currently only Human and Mouse supported"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Atropos_Trim", validator="Boolean", value="True"),
             my_controls.file_selector("--Anchored_Adapters_5p", tip="Full Path to Atropos 5' Adapter File"),
             my_controls.file_selector("--Anchored_Adapters_3p", tip="Full Path to Atropos 3' Adapter File"),
             my_controls.text_control("--Atropos_Aligner"),
             my_controls.restricted_selector("--NextSeq_Trim", validator="Boolean", value="True",
                                             tip="Specific Trimming for NextSeq."),
             my_controls.text_control("--Adapter_Mismatch_Fraction"),
             my_controls.int_control("--Read_Queue_Size", validator="Number"),
             my_controls.int_control("--Result_Queue_Size", validator="Number"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Demultiplex", validator="Boolean", value="False"),
             my_controls.text_control("--N_Limit", tip="Fraction of Sequence Allowed to be N"),
             my_controls.int_control("--Minimum_Length", validator="Number", tip="Minimum Read Length"),
             my_controls.int_control("--Trim5", validator="Number", tip="Additional Trimming of 5' End of FASTQ"),
             my_controls.int_control("--Trim3", validator="Number", tip="Additional Trimming of 3' End of FASTQ"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Aligner"),
             my_controls.text_control("--Aligner_Options", tip="Additional options such as -d 4"),
             my_controls.restricted_selector("--BWA_Method", tip="Only mem and aln supported"),

             my_controls.add_line(),
             my_controls.int_control("--Minimum_Family_Size", validator="Number"),
             my_controls.int_control("--UMT_Distance_Threshold", validator="Number",
                                     tip="For Deduplication, Generally Leave at 1"),
             my_controls.text_control("--Consensus_Freq_Threshold"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Target", tip="Genes or Exons"),
             my_controls.restricted_selector("--Filter_VCF_File", validator="Boolean", value="False"),
             my_controls.text_control("--SNP_Error_Freq", tip="For Filtering.  SNP Frequency in Error Matrix"),
             my_controls.int_control("--Boundary_Padding", validator="Number",
                                     tip="Number of Additional Nucleotides Outside of Regions"),
             my_controls.restricted_selector("--Strict_Boundaries", validator="Boolean", value="True"),
             my_controls.text_control("--Minimum_Allele_Freq",
                                      tip="When Filtering; What is the Minimum Minimum ALT Frequency Allowed"),
             my_controls.int_control("--Min_Fold_Increase", validator="Number"),
             my_controls.int_control("--Minimum_Allele_Count", validator="Number", tip="Minimum Alt Allele Count"),
             my_controls.int_control("--Minimum_Read_Depth", validator="Number"),
             my_controls.int_control("--Minimum_Base_Quality", validator="Number")
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class HaloPlexPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(HaloPlexPanel, self).__init__(parent, wx.ID_ANY, name="HaloPlexPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        self.panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        # title_icon = wx.StaticBitmap(self, wx.ID_DEFAULT, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Odin HaloPLEX Pipeline")
        title.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.FONTWEIGHT_BOLD))
        title.SetMinSize((int(parent.GetSize()[0] * 0.5), 45))
        title.SetForegroundColour("blue")
        title_sizer.Add(title, 0, wx.ALL)
        panel_sizer.Add(title_sizer, 0, wx.ALL | wx.CENTER)
        panel_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list = \
            [my_controls.file_selector("--Odin_API", tip="Full Path to Odin API"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.file_selector("--FASTQ1", tip="Full Path to FASTQ1"),
             my_controls.file_selector("--FASTQ2", tip="Full Path to FASTQ2"),
             my_controls.file_selector("--FASTQ3", tip="Full Path to FASTQ3"),
             my_controls.file_selector("--Ref_Seq", tip="Full Path to Genomic Reference Sequence"),
             my_controls.file_selector("--Aligner_Ref_Seq", tip="Full Path to Genomic Aligner Reference Sequence"),
             my_controls.file_selector("--Index_File"),
             my_controls.file_selector("--dbSNP", tip="Full Path to dbSNP file"),
             my_controls.file_selector("--Target_File", tip="Full Path to Genome Target File"),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Verbose", value="INFO"),
             my_controls.text_control("--Job_Name"),
             my_controls.int_control("--Spawn", validator="Number"),
             my_controls.restricted_selector("--Compression_Level", "9"),
             my_controls.restricted_selector("--Species", validator="Boolean",
                                             tip="Currently only Human and Mouse supported"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Atropos_Trim", validator="Boolean", value="True"),
             my_controls.file_selector("--Anchored_Adapters_5p", tip="Full Path to Atropos 5' Adapter File"),
             my_controls.file_selector("--Anchored_Adapters_3p", tip="Full Path to Atropos 3' Adapter File"),
             my_controls.text_control("--Atropos_Aligner"),
             my_controls.restricted_selector("--NextSeq_Trim", validator="Boolean", value="True",
                                             tip="Specific Trimming for NextSeq data."),
             my_controls.text_control("--Adapter_Mismatch_Fraction"),
             my_controls.int_control("--Read_Queue_Size", validator="Number"),
             my_controls.int_control("--Result_Queue_Size", validator="Number"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Demultiplex", validator="Boolean", value="False"),
             my_controls.text_control("--N_Limit", tip="Fraction of Sequence Allowed to be N"),
             my_controls.int_control("--Minimum_Length", validator="Number", tip="Minimum Read Length"),
             my_controls.int_control("--Trim5", validator="Number", tip="Additional Trimming of 5' End of FASTQ"),
             my_controls.int_control("--Trim3", validator="Number", tip="Additional Trimming of 3' End of FASTQ"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Aligner"),
             my_controls.text_control("--Aligner_Options", tip="Additional options such as -d 4"),
             my_controls.restricted_selector("--BWA_Method", tip="Only mem and aln supported"),

             my_controls.add_line(),
             my_controls.int_control("--Minimum_Family_Size", validator="Number"),
             my_controls.int_control("--UMT_Distance_Threshold", validator="Number",
                                     tip="For Deduplication, Generally Leave at 1"),
             my_controls.text_control("--Consensus_Freq_Threshold"),

             my_controls.add_line(),
             my_controls.restricted_selector("--Target", tip="Genes or Exons"),
             my_controls.restricted_selector("--Filter_VCF_File", validator="Boolean", value="False"),
             my_controls.text_control("--SNP_Error_Freq", tip="For Filtering.  SNP Frequency in Error Matrix"),
             my_controls.int_control("--Boundary_Padding", tip="Number of Additional Nucleotides Outside of Regions"),
             my_controls.restricted_selector("--Strict_Boundaries", validator="Boolean", value="True"),
             my_controls.text_control("--Minimum_Allele_Freq",
                                      tip="When Filtering; What is the Minimum Minimum ALT Frequency Allowed"),
             my_controls.int_control("--Min_Fold_Increase", validator="Number"),
             my_controls.int_control("--Minimum_Allele_Count", validator="Number", tip="Minimum Alt Allele Count"),
             my_controls.int_control("--Minimum_Read_Depth", validator="Number"),
             my_controls.int_control("--Minimum_Base_Quality", validator="Number")
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class WelcomePanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(WelcomePanel, self).__init__(parent, wx.ID_ANY, name="WelcomePanel")
        self.SetBackgroundStyle(wx.BG_STYLE_ERASE)
        self.box = wx.BoxSizer(wx.VERTICAL)

        m_text = wx.StaticText(self, wx.ID_ANY, "RTP Genomics L.L.C.")
        m_text.SetFont(wx.Font(20, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.BOLD))
        m_text.SetForegroundColour("blue")
        m_text.SetSize(m_text.GetBestSize())
        self.box.Add(m_text, 0, wx.ALL, 10)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.add_background)
        self.splash_image = parent.splash_image
        self.SetSizerAndFit(self.box)

    def add_background(self, event):
        """
         Add a picture to the background
         :param event:
         """
        # If the background file is missing the installation is likely damaged.
        if not os.path.isfile(self.splash_image):
            print("Background Image File {} not Found.  Installation Not Correct".format(self.splash_image))
            raise SystemExit(1)

        dc = wx.ClientDC(self)
        rect = self.GetUpdateRegion().GetBox()
        dc.SetClippingRegion(rect)

        bmp = wx.Bitmap(self.splash_image)
        screen_width, screen_height = self.GetClientSize()

        image_width = bmp.GetWidth()
        image_height = bmp.GetHeight()
        img = bmp.ConvertToImage()
        image_width = int((screen_width / image_width) * image_width)
        image_height = int((screen_height / image_height) * image_height)

        bmp = wx.Bitmap(img.Scale(image_width, image_height))

        x_pos = (screen_width - image_width) / 2
        y_pos = (screen_height - image_height) / 2
        dc.DrawBitmap(bmp, x_pos, y_pos)


def resource_path(relative_path):
    """
    Finds the splash screen image
    :param relative_path:
    :return:
    """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


class MainFrame(wx.Frame):
    """
    Sets up everything and displays our Welcome Screen
    """
    def __init__(self):

        # Setup a frame that is centered on the screen and opens at 75% of full screen.
        window_name = "Odin Bioinformatics GUI  v{}".format(__version__)
        display_x, display_y = wx.DisplaySize()
        scale = 0.75

        wx.Frame.__init__(self, None, wx.ID_ANY, title=window_name, size=(display_x * scale, display_y * scale))
        self.Centre()
        self.build_menu_bar()
        self.splash_image = resource_path("Splash_Image.jpg")
        self.pickle_file = "{0}{1}pickles{1}ODIN_parameters.pkl".format(Path(__file__).resolve().parent, os.sep)
        self.dataframe = self.dataframe_build()
        self.default_data = self.default_dict_build()
        print(self.GetSize(), wx.DisplaySize(), wx.version())

        self.welcome_panel = WelcomePanel(self)
        self.haloplex_panel = HaloPlexPanel(self)
        self.thruplex_panel = ThruPlexPanel(self)
        self.variant_search_panel = SkadiPanelVariantSearchPanel(self)
        self.filter_vcf_panel = SkadiFilterVCFPanel(self)
        self.build_error_matrix_panel = SkadiErrorMatrixPanel(self)
        self.annotate_vcf_panel = AnnotateVCFPanel(self)

        self.haloplex_panel.Hide()
        self.thruplex_panel.Hide()
        self.variant_search_panel.Hide()
        self.filter_vcf_panel.Hide()
        self.build_error_matrix_panel.Hide()
        self.annotate_vcf_panel.Hide()

        self.panel_dict = \
            {"HaloPlexPanel": self.haloplex_panel,
             "ThruPlexPanel": self.thruplex_panel,
             "WelcomePanel": self.welcome_panel,
             "SkadiPanelVariantSearchPanel": self.variant_search_panel,
             "SkadiFilterVCFPanel": self.filter_vcf_panel,
             "SkadiErrorMatrixPanel": self.build_error_matrix_panel,
             "AnnotateVCFPanel": self.annotate_vcf_panel
             }

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = self.welcome_panel
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)

        self.SetSizer(self.main_sizer)

    def dataframe_build(self):
        """

        :return:
        """
        try:
            with open(self.pickle_file, 'rb') as file:
                dataframe_dict = dill.load(file)

        except FileNotFoundError:
            # This might be the first time a user has run the program or the program has been reset.
            column_names = \
                ["--FASTQ1", "--FASTQ2", "--FASTQ3", "--Ref_Seq", "--Index_File", "--dbSNP", "--Target_File",
                 "--Working_Folder", "--Verbose", "--Job_Name", "--Spawn", "--Compression_Level", "--Species",
                 "--BAM_File", "--Atropos_Trim", "--Atropos_Aligner", "--Anchored_Adapters_5p", "--Anchored_Adapters_3p",
                 "--NextSeq_Trim", "--Adapter_Mismatch_Fraction", "--Read_Queue_Size", "--Result_Queue_Size",
                 "--Demultiplex", "--N_Limit", "--Minimum_Length", "--Trim5", "--Trim3", "--Aligner",
                 "--BWA_Method", "--Aligner_Options", "--Minimum_Family_Size", "--UMT_Distance_Threshold",
                 "--Consensus_Freq_Threshold", "--Target", "--Filter_VCF_File", "--Create_Control_Matrix",
                 "--Boundry_Padding", "--Strict_Boundaries", "--Min_Fold_Increase", "--Minimum_Allele_Count",
                 "--Minimum_Read_Depth", "--Minimum_Base_Quality", "--Atropos_Aligner", "--Aligner_Ref_Seq",
                 "--Variant_Search", "--SNP_Error_Freq", "--snpEff_Folder", "--Input_VCF", "--Java_Mem", "--Dataset",
                 "--Minimum_Allele_Freq", "--Odin_API", "--Options_File"
                 ]

            dataframe_dict = collections.defaultdict(list)
            for key in column_names:
                dataframe_dict[key] = []
        return dataframe_dict

    @staticmethod
    def default_dict_build():
        """
        This builds a dictionary with values that are library type and version specific.
        :return:
        """
        version1_default_dict = \
            {"--Compression_Level": ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
             "--Atropos_Trim": ["True", "False"], "--NextSeq_Trim": ["True", "False"],
             "--Read_Queue_Size": 500000, "--Demultiplex": ["True", "False"],
             "--Result_Queue_Size": 100000, "--BWA_Method": ["mem", "aln"], "--Minimum_Length": 25,
             "--Min_Fold_Increase": 3, "--Java_Mem": 12, "--Dataset": "GRCh38.86", "--Minimum_Allele_Freq": "0.001",
             "--Species": ["Human", "Mouse"], "--Variant_Target_Type": ["Exons", "Genes"], "--Aligner_Options": " ",
             "--Verbose": ["INFO", "DEBUG", "ERROR"], "--Spawn": 1, "--Filter_VCF_File": ["True", "False"],
             "--Create_Control_Matrix": ["True", "False"], "--Boundary_Padding": 20, "--N_Limit": "0.1",
             "--Strict_Boundaries": ["True", "False"], "--Minimum_Alt_Allele_Count": 1, "--Minimum_Read_Depth": 500,
             "--Minimum_Base_Quality": 20, "--Adapter_Mismatch_Fraction": "0.15", "--Consensus_Freq_Threshold": "0.55",
             "--UMT_Distance_Threshold": 1, "--Trim5": 0, "--Trim3": 0, '--Aligner': ["BWA", "BowTie2"],
             '--Minimum_Family_Size': 1, "--Minimum_Allele_Count": 1, "--Atropos_Aligner": "adapter",
             "--Target": ["Genes", "Exons"], "--Variant_Search": ["True", "False"], "--SNP_Error_Freq": "0.23",
             "--Include_All": ["True", "False"]
             }
        return version1_default_dict

    def build_menu_bar(self):
        """
        Build the menu bar
        """
        menu_bar = wx.MenuBar()

        # Build the "File" menu.
        file_menu = wx.Menu()
        welcome_app = file_menu.Append(wx.ID_ANY, "Home", "Go to Welcome Screen")
        save_app = file_menu.Append(wx.ID_SAVE, "&Save\tCtrl-S", "Save Parameter File")
        # run_app = file_menu.Append(wx.ID_EXECUTE, "&Run\tCtrl-R", "Save Parameter File and Run Program")
        file_menu.AppendSeparator()
        exit_app = file_menu.Append(wx.ID_EXIT, "E&xit\tAlt-x", "Close window and exit program.")
        menu_bar.Append(file_menu, "&File")

        # Build the "Tools" Menu
        tools_menu = wx.Menu()

        haloplex_app = tools_menu.Append(wx.ID_ANY, "Odin HaloPLEX")
        thruplex_app = tools_menu.Append(wx.ID_ANY, "Odin ThruPLEX")

        # Skadi SubMenu Items
        skadi_submenu = wx.Menu()
        variant_search_app = skadi_submenu.Append(wx.ID_ANY, "Variant_Search")
        filter_vcf_app = skadi_submenu.Append(wx.ID_ANY, "Filter_VCF (No Search)")
        build_error_matrix_app = skadi_submenu.Append(wx.ID_ANY, "Build Error Matrix")
        tools_menu.AppendSubMenu(skadi_submenu, "Skadi", )

        annotate_vcf_app = tools_menu.Append(wx.ID_ANY, "Annotate VCF")
        menu_bar.Append(tools_menu, "&Tools")

        # Build the "Help" menu
        help_menu = wx.Menu()
        about_item = help_menu.Append(wx.ID_ABOUT)
        menu_bar.Append(help_menu, "&Help")
        self.SetMenuBar(menu_bar)

        # Bind Menu Items to Actions
        self.Bind(wx.EVT_MENU, self.save_parameter_file, save_app)
        # self.Bind(wx.EVT_MENU, self.on_run, run_app)
        self.Bind(wx.EVT_MENU, self.on_exit, exit_app)
        self.Bind(wx.EVT_MENU, self.on_about, about_item)

        # Panel switching
        self.Bind(wx.EVT_MENU, lambda event, temp=haloplex_app:
                  self.switch_panels(event, self.haloplex_panel.GetName()), haloplex_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=thruplex_app:
                  self.switch_panels(event, self.thruplex_panel.GetName()), thruplex_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=variant_search_app:
                  self.switch_panels(event, self.variant_search_panel.GetName()), variant_search_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=filter_vcf_app:
                  self.switch_panels(event, self.filter_vcf_panel.GetName()), filter_vcf_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=build_error_matrix_app:
                  self.switch_panels(event, self.build_error_matrix_panel.GetName()), build_error_matrix_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=welcome_app:
                  self.switch_panels(event, self.welcome_panel.GetName()), welcome_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=annotate_vcf_app:
                  self.switch_panels(event, self.annotate_vcf_panel.GetName()), annotate_vcf_app)

    def switch_panels(self, event, switch_id):
        """
        Handles changing windows in GUI
        :param event:
        :param switch_id:
        """
        self.main_sizer.Detach(self.panel)
        self.panel.Hide()
        self.panel = self.panel_dict[switch_id]
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)
        self.panel.Show()
        self.panel.Fit()
        self.Layout()

    def on_exit(self, event):
        """Close the frame, terminating the application.
        :param event:
        """
        self.Close(True)

    def on_about(self, event):
        """
        Creates the About window
        :param event:
        """
        about_info = wx.adv.AboutDialogInfo()

        about_info.SetName("Odin Bioinformatics App")
        about_info.SetCopyright(__copyright__)

        about_info.SetDescription(wordwrap(
            "GUI for Odin.\nVersion: {}".format(__version__), 350, wx.ClientDC(self)))
        about_info.SetDevelopers([__author__])
        # about_info.License = wordwrap("Completely and totally open source!", 500, wx.ClientDC(self))
        wx.adv.AboutBox(about_info)

    def on_run(self, event):
        """
        Currently not used.
        :param event:
        """
        shellfile_name = self.save_parameter_file("Run")
        print(shellfile_name)
        if shellfile_name is not None:
            subprocess.run([shellfile_name], shell=True)

    def save_parameter_file(self, event):
        """

        :param event:
        :return:
        """
        def outstring_build():
            """
            :return:
            """
            file_body = "--HaloPLEX\t{}\n" \
                        "--ThruPLEX\t{}\n" \
                        "--Skadi\t{}\n" \
                        "--Annotate_Results\t{}\n" \
                        "--Create_Error_Matrix\t{}\n" \
                        "--FASTQ_Quality\tFalse\n" \
                        "{}\n".format(haloplex, thruplex, skadi, annotate, error_matrix, additional_parameters)
            working_folder = ""
            job = ""
            api = ""
            options = ""

            for i in range(len(panel.control_dict)):
                d = panel.control_dict[i][1]
                # Validate data.  If validation fails this will take us back to the form.
                with suppress(AttributeError):
                    if not d.GetValidator().Validate(d):
                        return job, working_folder, file_body, False

                ctrl_name = d.GetName()
                if ctrl_name == 'line':
                    file_body += "\n"
                    continue
                ctrl_value = str(d.GetValue()).strip()

                # Linux does not add the trailing slash to folders.  I find that confusing.
                if ctrl_name == "--Odin_API":
                    api = ctrl_value
                elif ctrl_name == "--Working_Folder":
                    working_folder = ctrl_value
                    ctrl_value += "/"
                elif ctrl_name == "--Job_Name":
                    job = ctrl_value
                elif ctrl_name == "--NextSeq_Trim":
                    if ctrl_value == "True":
                        ctrl_value = 1
                    else:
                        ctrl_value = 0
                elif ctrl_name == "--Options_File":
                    options = ctrl_value

                if not ctrl_name == "--Odin_API" or not ctrl_name == "--Options_File":
                    file_body += "{}\t{}\n".format(ctrl_name, ctrl_value)
                self.dataframe[ctrl_name].append(d.GetValue())
                self.dataframe[ctrl_name] = list(set(self.dataframe[ctrl_name]))

            return job, file_body, True, api, options

        haloplex = "False"
        thruplex = "False"
        skadi = "False"
        annotate = "False"
        error_matrix = "False"
        submodule = "Configuration Error"
        additional_parameters = "\n"

        if self.haloplex_panel.IsShown():
            submodule = "HaloPLEX Pipeline"
            haloplex = "True"
            panel = self.haloplex_panel

        elif self.thruplex_panel.IsShown():
            submodule = "ThruPLEX Pipeline"
            thruplex = "True"
            panel = self.thruplex_panel

        elif self.variant_search_panel.IsShown():
            submodule = "Skadi Variant Search"
            skadi = "True"
            panel = self.variant_search_panel
            additional_parameters += "--Variant_Search\tTrue\n"

        elif self.filter_vcf_panel.IsShown():
            submodule = "Skadi Filter VCF File"
            skadi = "True"
            panel = self.filter_vcf_panel
            additional_parameters += "--Variant_Search\tFalse\n" \
                                     "--Filter_VCF_File\tTrue\n"

        elif self.build_error_matrix_panel.IsShown():
            error_matrix = "True"
            submodule = "Skadi Build Error Matrix"
            panel = self.build_error_matrix_panel

        elif self.annotate_vcf_panel.IsShown():
            annotate = "True"
            submodule = "Annotate VCF File"
            panel = self.annotate_vcf_panel

        job_name, parameter_body, validation_pass, api_path, options_path = outstring_build()

        date = datetime.datetime.now().strftime("%Y%m%d")
        shebang = "#!/bin/bash\n" \
                  "#Parameter file to run Odin Pipeline module {}\n" \
                  "#File generated {}\n\n".format(submodule, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        outfile_name = "{}{}run_{}_{}.sh".format(options_path, os.sep, job_name, date)
        cmd = "python3 {} --options_file {}\nexit\n\n".format(api_path, outfile_name)

        # Save the current choices in the pickles folder.
        with open(self.pickle_file, 'wb') as file:
            dill.dump(self.dataframe, file, protocol=-1)

        # Write the parameter file.
        outstring = "{}{}{}".format(shebang, cmd, parameter_body)
        filename = "{}_{}.sh".format(job_name, date)

        if event == "Run" and validation_pass:
            parameter_file = open(outfile_name, 'w')
            parameter_file.write(outstring)
            parameter_file.close()
            return outfile_name

        with wx.FileDialog(self, "Save Odin Options File",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT, defaultFile=outfile_name) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # the user changed their mind

            # save the current contents in the file
            pathname = fileDialog.GetPath()

            try:
                parameter_file = open("{}".format(pathname), 'w')
                parameter_file.write(outstring)
                parameter_file.close()
            except IOError:
                wx.LogError("Cannot save {} in {}".format(filename, pathname))


def main():
    app = wx.App()
    frame = MainFrame()
    frame.Show()
    app.MainLoop()
    # wx.lib.inspection.InspectionTool().Show()


if __name__ == '__main__':
    main()
