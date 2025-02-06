
import React, { useRef, useState } from "react";
import { connect } from "react-redux";

import { Modal } from "antd";
import Draggable from "react-draggable";
import cls from "./index.css";

const defaultBounds = { left: 0, top: 0, bottom: 0, right: 0 };

const ConversationModal = () => {
  const draggleRef = useRef(null);
  const [bounds, setBounds] = useState(defaultBounds);
  const [visible, setVisible] = useState(true);

  const closeModal = () => setVisible(false);

  const onStart = (_event, uiData) => {
    const { clientWidth, clientHeight } = window.document.documentElement;
    const targetRect = draggleRef.current?.getBoundingClientRect();
    if (!targetRect) {
      return;
    }
    setBounds({
      left: -targetRect.left + uiData.x,
      right: clientWidth - (targetRect.right - uiData.x),
      top: -targetRect.top + uiData.y,
      bottom: clientHeight - (targetRect.bottom - uiData.y),
    });
  };

  const modalRender = (modal) => (
    <Draggable
      bounds={bounds}
      nodeRef={draggleRef}
      onStart={onStart}
      defaultClassName={cls.mask}
    >
      <div ref={draggleRef}>{modal}</div>
    </Draggable>
  );

  return (
    <Modal
      open={visible}
      onCancel={closeModal}
      mask={false}
      maskClosable={false}
      classNames={{ mask: cls.mask, wrapper: cls.mask }}
      modalRender={modalRender}
      footer={null}
      title="RegionName123 (umap)567cells"
    >
      <div>Hello World</div>
    </Modal>
  );
};

const mapStateToProps = (state) => ({
  selections: state.graphSelection.selections,
});

export default connect(mapStateToProps)(ConversationModal);
