import React, {PureComponent} from 'react';

export class Popup extends PureComponent {
  render() {
    const { title, x, y } = this.props;
    const style = {
      position: 'absolute',
      top: y,
      left: x,
      maxWidth: '200px',
      padding: '10px',
      color: 'white',
      backgroundColor: 'black',
      pointerEvents: 'none',
      transform: 'translate(10px, -50%)',
    };

    const arrowStyle = {
      position: 'absolute',
      top: '50%',
      left: '-14px',
      width: '7px',
      height: '5px',
      boxSizing: 'border-box',
      transform: 'translateY(-50%)',
      border: '7px solid transparent',
      borderRight: '7px solid black',
    };

    return (
      <div style={style}>
        <div style={arrowStyle} />
        {title}
      </div>
    )
  };
}

Popup.propTypes = {
  id: React.PropTypes.string,
  title: React.PropTypes.string,
  x: React.PropTypes.number,
  y: React.PropTypes.number,
};

Popup.defaultProps = {
  id: 'id',
  title: '',
  x: 0,
  y: 0,
};
