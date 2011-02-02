#! /usr/local/bin/python

# communicate.py
# Demo widget code from http://zetcode.com/wxpython/firststeps/
# To run type python2.7 communicate.py &

import wx
import Image


class LeftPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id, style=wx.BORDER_SUNKEN)

        self.text = parent.GetParent().rightPanel.text

        button_run = wx.Button(self, -1, 'Run', (10, 10))
        button_next = wx.Button(self, -1, 'Next', (10, 60))
        button_back = wx.Button(self, -1, 'Back', (10, 60))

        self.Bind(wx.EVT_BUTTON, self., id=button_run.GetId())
        self.Bind(wx.EVT_BUTTON, self.OnPlus, id=button_next.GetId())
        self.Bind(wx.EVT_BUTTON, self.OnMinus, id=button_back.GetId())

    def OnPlus(self, event):
        value = int(self.text.GetLabel())
        value = value + 1
        self.text.SetLabel(str(value))

    def OnMinus(self, event):
        value = int(self.text.GetLabel())
        value = value - 1
        self.text.SetLabel(str(value))


class RightPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id, style=wx.BORDER_SUNKEN)
        self.text = wx.StaticText(self, -1, '0', (40, 60))


class Communicate(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(280, 200))

        panel = wx.Panel(self, -1)
        self.rightPanel = RightPanel(panel, -1)

        leftPanel = LeftPanel(panel, -1)

        hbox = wx.BoxSizer()
        hbox.Add(leftPanel, 1, wx.EXPAND | wx.ALL, 5)
        hbox.Add(self.rightPanel, 1, wx.EXPAND | wx.ALL, 5)

        panel.SetSizer(hbox) 
        self.Centre()
        self.Show(True)

app = wx.App()
Communicate(None, -1, 'widgets communicate')
app.MainLoop()

